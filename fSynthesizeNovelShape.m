%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jun/2017 - XXX/2018
%Funcion: fSynthesizeNovelShape
%Input:
%   segmentsArray:      structure array containing the segments that were previously computed for all 
%                       the input shapes
%   samplingsArray:     structure array containing the samples of the original Shapes
%   templateIndices:    indices of the segments of the current Template shape in the array 
%                       "segmentsArray", (see fGetSegIndicesPerShape)
%   selSegs:            indices of the segments that were selected as replacements for each segment 
%                       of the template, (see fSampleSegsForNewShape)
%   segmentsAdj:        structure containing graphs that represent the connectivity 
%                       of the segments for each Shape (see fGetAdjacencyPerShape)
%   pathToSave:         full path (including file name and extension) for
%                       saving the new shape as a point cloud, in a PLY file
%   doLocAlign:         (OPTIONAL) indicates if we apply or not the process of local
%                       alignment between adjacent shapes
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%   newSegments:    (OPTIONAL) array of structures containing the info of the new
%                   segments used for the synthesized shape, i.e, the points and the normals
%   transforms:     array of structures containing the transformation that was
%                   applied to each new (replacement) segment. This array contains 
%                   a rotation matrix, a scaling matrix, and a translation vertex 
%                   per segment
%
% Function that synthesizes a new shape (a point cloud) using the new segments ("selSegs") to replace 
% each segment of a Template Shape, and saves the result to a PLY file in
% the folder given by "pathToSave"
%
% (See also fSampleSegsForNewShape, fGetAdjacencyPerShape, fSavePly_FromSegments)
%%%%%
%(1) References: 
%   (a) Zheng, Y., et al. "Smart Variations: Functional Substructures for Part Compatibility" (2013)
%   (b) Prakhar, J., et al. "Assembly-based conceptual 3D modeling with unlabeled components using probabilistic factor graph" (2016)
%   (c) Cormen TH. "Introduction to Algorithms" (2009)
%   (d) https://en.wikipedia.org/wiki/Breadth-first_search
%%%%%
%--------------------------------------------------------------------------
function [err, newSegments, transforms] = fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, selSegs, segmentsAdj, pathToSave, doLocAlign)
    try
        if length(selSegs) ~= length(templateIndices)
            sErr= 'Error: the number of segments of the template must match the number of the new segments. Check the input parameters';
            errordlg(sErr);
            return;
        end
        if (~exist('doLocAlign', 'var'))
            doLocAlign= true;
        end
        
        newSegments= struct([]);
        transforms= struct([]);
        colorPerSeg= zeros(length(templateIndices), 3);
        templateSegDone= false(length(templateIndices), 1);
        neighIsSymm= false(length(templateIndices), 1);
        
        %Initialize the array of structures of new segments
        for i= 1 : length(templateIndices)
            nS.points = [];
            nS.normals= [];
            newSegments= [newSegments nS];
            
            tra.rot= zeros(3);
            tra.scale= zeros(3);
            tra.translat= zeros(1, 3);
            tra.reflec= false;
            transforms= [transforms tra];
        end
        
        %I. Replace the segments of the Template, trying to match the size
        %and pos. of each segment, to create an initial novel Shape...
        shapeIxR= segmentsArray(templateIndices(1)).ixShape;
        centShapeRef= mean(samplingsArray(shapeIxR).samples);
        for i= 1 : size(templateIndices, 2)            
            if (~templateSegDone(i))
                swc= false;
                doCalcRot= true;                
                iSel= i;
                
                newSeg= segmentsArray(selSegs(iSel));
                ptsTemp= samplingsArray(newSeg.ixShape).samples(newSeg.ixSamples, :);
                normsNew= samplingsArray(newSeg.ixShape).normals(newSeg.ixSamples, :);
                
                refSeg= segmentsArray(templateIndices(iSel));
                ptsRef= samplingsArray(shapeIxR).samples(refSeg.ixSamples, :);
                
                %%%%%
                %For the case when we preserve some segments of the original template, no transformation is required
                if (templateIndices(iSel) == selSegs(iSel))
                    newSegments(iSel).points= ptsTemp;
                    newSegments(iSel).normals= normsNew;                    
                    
                    transforms(iSel).rot= eye(3);
                    transforms(iSel).scale= eye(3);
                    
                    colorPerSeg(iSel, :)= newSeg.color;
                    templateSegDone(iSel)= true;
                    continue;
                end
                %%%%%
                
                %Find the label of the symmetric segments w.r.t to the current template's segment
                [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(i), segmentsArray, templateIndices);
                if (err)
                    break;
                end
                                
                %Find the label of the symmetric segments w.r.t to the new (replacement) segment
                [err, symGroupNew, ixSymGN] = fGetSymmetryChain(selSegs(i), segmentsArray);
                if (err)
                    break;
                end                

                if (~isempty(symGroupNew))
%                     if (length(symGroupTemplate) == length(symGroupNew)) && ((length(symGroupTemplate) == 2) || (length(symGroupTemplate) == 4))
                    if ((length(symGroupNew) == 2) || (length(symGroupNew) == 4)) && ((length(symGroupTemplate) == 2) || (length(symGroupTemplate) == 4))
                        centRef= mean(ptsRef);
                        centNew= mean(ptsTemp);
                        centShapeNew= mean(samplingsArray(newSeg.ixShape).samples);
                        
                        if (length(symGroupTemplate) == 2)
                            if ((centRef(3) < centShapeRef(3)) && (centNew(3) < centShapeNew(3))) || ...
                                ((centRef(3) >= centShapeRef(3)) && (centNew(3) >= centShapeNew(3)))
                                iSel= i;                                
                            else
                                swc= true;
                                iSel= symGroupTemplate(2);
                                refSeg= segmentsArray(templateIndices(iSel));
                                ptsRef= samplingsArray(shapeIxR).samples(refSeg.ixSamples, :);
                            end
                        elseif (length(symGroupTemplate) == 4)
                            for iisl= 1 : 4
                                refSegCan= segmentsArray(templateIndices(symGroupTemplate(iisl)));
                                ptsRefCan= samplingsArray(shapeIxR).samples(refSegCan.ixSamples, :);
                                centRefCan= mean(ptsRefCan);                                
                                if (all(sign(centRefCan) == sign(centNew)))
                                    iSel= symGroupTemplate(iisl);
                                    swc= (iSel ~= i);
                                    refSeg= refSegCan;
                                    ptsRef= ptsRefCan;
                                    break;
                                end
                            end
                        end
                    end                    
                end
                
                szRef= max(ptsRef)-min(ptsRef);
                obbR= part_obb(ptsRef);
                
                obbN= part_obb(ptsTemp);
                
                %Find the neighbours of the template's segment
                ixNeighsTemplate= find(segmentsAdj(shapeIxR).adjacencyMat(iSel, :) == 1);
                
                %Find the neighbours of the replacement segment
                ixNeighsNew= find(segmentsAdj(newSeg.ixShape).adjacencyMat(newSeg.segLabel, :) == 1);
                
                %%%
                %The following lines of code are for a very particular case: when a pair of segments are neighbours,
                %and also are symmetric, we can determine the rotation in a more direct way
                
                %Check if any symmetric segment is also a neighbour
                ixNeighAlsoSymR= 0;
                ixNeighAlsoSymN= 0;
                refSymIsNeigh= find(ismember(symGroupTemplate, ixNeighsTemplate));
                newSymIsNeigh= find(ismember(symGroupNew, ixNeighsNew));
                
                if (~isempty(refSymIsNeigh))
                    if (~isempty(newSymIsNeigh))                        
                        %%%
                        for nsym= 1 : length(refSymIsNeigh)
                            if (segmentsArray(templateIndices(iSel)).segmentMainAx == ...
                                    segmentsArray(templateIndices(symGroupTemplate(refSymIsNeigh(nsym)))).segmentMainAx)
                                ixNeighAlsoSymR= nsym;
                                break;
                            end
                        end
                        
                        for nsym= 1 : length(newSymIsNeigh)                            
%                             if (newSeg.segmentMainAx == segmentsArray(symGroupNew(newSymIsNeigh(nsym))).segmentMainAx)
                            if (newSeg.segmentMainAx == segmentsArray(ixSymGN(newSymIsNeigh(nsym))).segmentMainAx)
                                ixNeighAlsoSymN= nsym;
                                break;
                            end
                        end
                        %%%
                        
                        if (ixNeighAlsoSymR > 0) && (ixNeighAlsoSymN > 0)
                            if (segmentsArray(templateIndices(iSel)).segmentMainAx == ...
                                    segmentsArray(templateIndices(symGroupTemplate(refSymIsNeigh(ixNeighAlsoSymR)))).segmentMainAx) && ...
                                    (newSeg.segmentMainAx == segmentsArray(ixSymGN(newSymIsNeigh(ixNeighAlsoSymN))).segmentMainAx) && ...
                                    (segmentsArray(templateIndices(iSel)).segmentMainAx == newSeg.segmentMainAx)
                                %%%
                                %For these special segments, we save a flag that is used
                                %later, when aligning the new neighbor segments
                                neighIsSymm(iSel)= true;
                                neighIsSymm(symGroupTemplate(refSymIsNeigh))= true;
                                %%%

                                doCalcRot= false;
                                centRef= mean(ptsRef);
                                centNew= mean(ptsTemp);
                                centShapeNew= mean(samplingsArray(newSeg.ixShape).samples);

                                if (segmentsAdj(shapeIxR).adjacencyInfo(iSel, symGroupTemplate(refSymIsNeigh(ixNeighAlsoSymR))).side == ...
                                        segmentsAdj(newSeg.ixShape).adjacencyInfo(newSeg.segLabel, symGroupNew(newSymIsNeigh(ixNeighAlsoSymN))).side)
                                    %Same side
                                    if ((centRef(1) < centShapeRef(1)) && (centNew(1) < centShapeNew(1))) || ...
                                            ((centRef(1) >= centShapeRef(1)) && (centNew(1) >= centShapeNew(1)))
                                        mR= eye(3);
                                    elseif (centRef(1) < centShapeRef(1)) && (centNew(1) >= centShapeNew(1))
                                        %Rotate about global upright vector
                                        if (~swc)
                                            mR= vrrotvec2mat([[0 1 0] -pi/2]);
                                        else
                                            mR= vrrotvec2mat([[0 1 0] pi/2]);
                                        end
                                    elseif (centRef(1) >= centShapeRef(1)) && (centNew(1) < centShapeNew(1))
                                        %Rotate about global upright vector
                                        if (~swc)
                                            mR= vrrotvec2mat([[0 1 0] pi/2]);
                                        else
                                            mR= vrrotvec2mat([[0 1 0] -pi/2]);
                                        end
                                    end
                                else
                                    %Different side
                                    if ((centRef(1) < centShapeRef(1)) && (centNew(1) < centShapeNew(1))) || ...
                                            ((centRef(1) >= centShapeRef(1)) && (centNew(1) >= centShapeNew(1)))
                                        %Reflect points through xy plane
                                        mR= [1 0 0 ; 0 1 0 ; 0 0 -1];
                                        transforms(iSel).reflec= true;
                                    else
                                        %Rotate pi using global upright vector
                                        mR= vrrotvec2mat([[0 1 0] pi]);
                                    end
                                end
                            end
                        end
                    end
                end
                %%%
                
                mainAxOfRef= refSeg.segmentMainAx;
                mainAxNew= newSeg.segmentMainAx;
                if (doCalcRot)
                    axRot= cross(obbR.axes(1:3, mainAxOfRef), obbN.axes(1:3, mainAxNew));
                    dp= sum(obbR.axes(1:3, mainAxOfRef) .* obbN.axes(1:3, mainAxNew));
                    if (dp > 1)
                        dp= 1;
                    elseif (dp < -1)
                        dp= -1;
                    end
                    angRot= acos(dp);
                    if (dp < 0)
                        axRot= -axRot;
                    elseif (dp == 0)
                        aConsR= zeros(3, 1);
                        aConsR(mainAxOfRef)= 1;
                        aConsN= zeros(3, 1);
                        aConsN(mainAxNew)= 1;
                        if (aConsR == aConsN)
                            axRot= aConsN;
                        else
                            axRot= cross(aConsR, aConsN);
                        end 
                    end
                    mR= vrrotvec2mat([axRot' angRot]);
                    
                    ptsNew= ptsTemp * mR;
                    ptsNew= fTransAndScalePoints(ptsNew, ptsRef, szRef);
                    
                    %Compute a value similar to the RMSE between the two set of points 
                    [~, minD]= knnsearch(ptsNew, ptsRef);
                    minRmse= sqrt(sum(minD.^2)/length(minD));
                    
                    %We try several other possible rotations...
                    axRot= cross(obbR.axes(1:3, mainAxNew), obbN.axes(1:3, mainAxOfRef));
                    dp= sum(obbR.axes(1:3, mainAxNew) .* obbN.axes(1:3, mainAxOfRef));
                    if (dp > 1)
                        dp= 1;
                    elseif (dp < -1)
                        dp= -1;
                    end
                    angRot= acos(dp);
                    if (dp < 0)
                        axRot= -axRot;
                    elseif (dp == 0)
                        aConsR= zeros(3, 1);
                        aConsR(mainAxOfRef)= 1;
                        aConsN= zeros(3, 1);
                        aConsN(mainAxNew)= 1;
                        if (aConsR == aConsN)
                            axRot= aConsN;
                        else
                            axRot= cross(aConsR, aConsN);
                        end 
                    end
                    tmpMR= vrrotvec2mat([axRot' angRot]);

                    ptsNew= ptsTemp * tmpMR;
                    
                    ptsNew= fTransAndScalePoints(ptsNew, ptsRef, szRef);
                    
                    %Compute a value similar to the RMSE between the two set of points and
                    %compare it with the previous one
                    [~, minD]= knnsearch(ptsNew, ptsRef);
                    tmpRmse= sqrt(sum(minD.^2)/length(minD));
                    if (tmpRmse < minRmse)
                        minRmse= tmpRmse;
                        %We keep the rotation having the min. error
                        mR= tmpMR;
                    end
              
                    axRot= cross(obbR.axes(1:3, mainAxNew), obbN.axes(1:3, mainAxNew));
                    dp= sum(obbR.axes(1:3, mainAxNew) .* obbN.axes(1:3, mainAxNew));
                    if (dp > 1)
                        dp= 1;
                    elseif (dp < -1)
                        dp= -1;
                    end
                    angRot= acos(dp);
                    if (dp < 0)
                        axRot= -axRot;
                    elseif (dp == 0)
                        aConsR= zeros(3, 1);
                        aConsR(mainAxOfRef)= 1;
                        aConsN= zeros(3, 1);
                        aConsN(mainAxNew)= 1;
                        if (aConsR == aConsN)
                            axRot= aConsN;
                        else
                            axRot= cross(aConsR, aConsN);
                        end 
                    end
                    tmpMR= vrrotvec2mat([axRot' angRot]);

                    ptsNew= ptsTemp * tmpMR;
                    
                    ptsNew= fTransAndScalePoints(ptsNew, ptsRef, szRef);
                    
                    %Compute a value similar to the RMSE between the two set of points and
                    %compare it with the previous one
                    [~, minD]= knnsearch(ptsNew, ptsRef);
                    tmpRmse= sqrt(sum(minD.^2)/length(minD));
                    if (tmpRmse < minRmse)
                        minRmse= tmpRmse;
                        %We keep the rotation having the min. error
                        mR= tmpMR;
                    end

                    axRot= cross(obbR.axes(1:3, mainAxOfRef), obbN.axes(1:3, mainAxOfRef));
                    dp= sum(obbR.axes(1:3, mainAxOfRef) .* obbN.axes(1:3, mainAxOfRef));
                    if (dp > 1)
                        dp= 1;
                    elseif (dp < -1)
                        dp= -1;
                    end
                    angRot= acos(dp);
                    if (dp < 0)
                        axRot= -axRot;
                    elseif (dp == 0)
                        aConsR= zeros(3, 1);
                        aConsR(mainAxOfRef)= 1;
                        aConsN= zeros(3, 1);
                        aConsN(mainAxNew)= 1;
                        if (aConsR == aConsN)
                            axRot= aConsN;
                        else
                            axRot= cross(aConsR, aConsN);
                        end 
                    end
                    tmpMR= vrrotvec2mat([axRot' angRot]);

                    ptsNew= ptsTemp * tmpMR;
                    
                    ptsNew= fTransAndScalePoints(ptsNew, ptsRef, szRef);
                    
                    %Compute a value similar to the RMSE between the two set of points and
                    %compare it with the previous one
                    [~, minD]= knnsearch(ptsNew, ptsRef);
                    tmpRmse= sqrt(sum(minD.^2)/length(minD));
                    if (tmpRmse < minRmse)
                        %We keep the rotation having the min. error
                        mR= tmpMR;
                    end
                    
                    if (mainAxOfRef == mainAxNew)
                        %If segments have the same main axis, try also without doing any rotation
                        tmpMR= eye(3);

                        ptsNew= ptsTemp * tmpMR;

                        ptsNew= fTransAndScalePoints(ptsNew, ptsRef, szRef);

                        %Compute a value similar to the RMSE between the two set of points and
                        %compare it with the previous one
                        [~, minD]= knnsearch(ptsNew, ptsRef);
                        tmpRmse= sqrt(sum(minD.^2)/length(minD));
                        if (tmpRmse < minRmse)
                            %We keep the rotation having the min. error
                            mR= tmpMR;
                        end
                    end
                end
                
                %We apply the rotation to the points of the new segment
                ptsNew= ptsTemp * mR;
                
                %We scale the new segment according to the proportions of the template's segment
                [ptsNew, mScale]= fTransAndScalePoints(ptsNew, ptsRef, szRef);

                %We save the points of the new segment
                newSegments(iSel).points= ptsNew;
                newSegments(iSel).normals= normsNew * mR;
                newSegments(iSel).normals= bsxfun(@rdivide, newSegments(iSel).normals, sqrt(sum(newSegments(iSel).normals .^ 2, 2)));
                %We save also the rotation and scaling applied. (We'll get and save the translation at the end)
                transforms(iSel).rot= mR;
                transforms(iSel).scale= mScale;
                
                if (~isempty(symGroupTemplate))
                    %Now we look for the segments that are symmetric, and we reflect (rotate) 
                    %the segment that we already transformed
                    symGroup= symGroupTemplate(symGroupTemplate ~= iSel);
                    cenB= mean(ptsRef);
                    for sy= 1 : length(symGroup)
                        transforms(symGroup(sy)).scale= mScale;
                        
                        tMRef= eye(3);
                        ptsSym= samplingsArray(shapeIxR).samples(segmentsArray(templateIndices(symGroup(sy))).ixSamples, :);
                        cenSym= mean(ptsSym);
                        
                        %Find the symmetry (reflection) plane for the symmetric segment
                        tempV= (cenB-cenSym)';
                        tempV= tempV / norm(tempV);
                        maxDp= 0;
                        selAx= 0;
                        for a= 1 : 3
                           axR= zeros(3, 1);
                           axR(a)= 1;
                           dp= abs(sum(tempV .* axR));
                           if (dp > maxDp && dp > 0.95)
                              maxDp= dp;
                              selAx= a; 
                           end
                        end
                        if (selAx >0)
                            tMRef(:, selAx)= -tMRef(:, selAx);
                            transforms(symGroup(sy)).reflec= true;
                        else
                            %We try rotating the points by pi using the global upright vector
                            tMRef= vrrotvec2mat([[0 1 0] pi]);
                        end
                        
                        %Reflect the points of the replacement segment through the reflection plane
                        ptsNS= newSegments(iSel).points * tMRef;

                        if (~neighIsSymm(iSel))
                            deltaNS= cenSym - mean(ptsNS);
                            newSegments(symGroup(sy)).points= ptsNS + repmat(deltaNS, size(ptsNS, 1), 1);
                            delta2= mean(newSegments(iSel).points) - mean(newSegments(symGroup(sy)).points);
                            if (selAx >0)
                                delta2(selAx)= 0;
                            else
                                delta2(1)= 0;
                                delta2(3)= 0;
                            end
                            newSegments(symGroup(sy)).points= newSegments(symGroup(sy)).points + repmat(delta2, size(newSegments(symGroup(sy)).points, 1), 1);
                        else
                            if (selAx == 3)
                                deltaNS= zeros(1, 3);
                                if (sign(cenB(3)) == -1)
                                    difAxZ= abs(max(newSegments(iSel).points) - min(ptsNS));
                                    if (difAxZ(3) > 0.004)
                                        deltaNS(3)= (difAxZ(3) / 2.0);
                                        newSegments(iSel).points= newSegments(iSel).points + repmat(deltaNS, size(ptsNS, 1), 1);
                                        newSegments(symGroup(sy)).points= ptsNS + repmat(-deltaNS, size(ptsNS, 1), 1);
                                    else
                                        newSegments(symGroup(sy)).points= ptsNS;
                                    end
                                elseif (sign(cenB(3)) == 1)
                                    difAxZ= abs(min(newSegments(iSel).points) - max(ptsNS));
                                    if (difAxZ(3) > 0.004)
                                        deltaNS(3)= -(difAxZ(3) / 2.0);
                                        newSegments(iSel).points= newSegments(iSel).points + repmat(deltaNS, size(ptsNS, 1), 1);
                                        newSegments(symGroup(sy)).points= ptsNS + repmat(-deltaNS, size(ptsNS, 1), 1);
                                    else
                                        newSegments(symGroup(sy)).points= ptsNS;
                                    end
                                end
                            else
                                newSegments(symGroup(sy)).points= ptsNS;
                            end
                        end
                        %We save the normals of the new segment
                        newSegments(symGroup(sy)).normals= newSegments(iSel).normals * tMRef;
                        newSegments(symGroup(sy)).normals= bsxfun(@rdivide, newSegments(symGroup(sy)).normals, ...
                            sqrt(sum(newSegments(symGroup(sy)).normals .^ 2, 2)));
                        %We save the new rotation applied
                        transforms(symGroup(sy)).rot= mR * tMRef;
                    end
                    templateSegDone(symGroupTemplate)= true;
                    colorPerSeg(symGroupTemplate, :)= repmat(newSeg.color, length(symGroupTemplate), 1);
                else
                    templateSegDone(i)= true;
                    colorPerSeg(i, :)= newSeg.color;
                end                
            end
            
            if(~doLocAlign)
                %Calculate and save the translation that was applied
                refSeg= segmentsArray(templateIndices(i));
                ptsRef= samplingsArray(shapeIxR).samples(refSeg.ixSamples, :);

                newSegTr= segmentsArray(selSegs(i));
                ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                ptsNewTr= ptsNewTr * transforms(i).rot;
                ptsNewTr= ptsNewTr * transforms(i).scale;

                deltaP= mean(ptsRef) - mean(ptsNewTr);
                transforms(i).translat= deltaP;
            end
        end
        
        %%%
        %Fix for the orientation of some segments that are marked first
        %incorrectly for being "reflected"
        for i= 1 : length(templateIndices)
            [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(i), segmentsArray, templateIndices);
            if (err)
                break;
            end
            if (length(symGroupTemplate) ==  4) && (neighIsSymm(i))
                newSeg= segmentsArray(selSegs(i));
                ptsTemp= samplingsArray(newSeg.ixShape).samples(newSeg.ixSamples, :);
                centNew= mean(ptsTemp);                        
                for iisl= 1 : 4
                    refSegCan= segmentsArray(templateIndices(symGroupTemplate(iisl)));
                    ptsRefCan= samplingsArray(shapeIxR).samples(refSegCan.ixSamples, :);
                    centRefCan= mean(ptsRefCan);                                
%                     if (all(sign(centRefCan) == sign(centNew)))
                    if (sign(centRefCan(1)) == sign(centNew(1))) && (sign(centRefCan(3)) == sign(centNew(3)))
                        iSelChk= symGroupTemplate(iisl);
                        ptsRef= ptsRefCan;
                        break;
                    end
                end
                
                if (transforms(iSelChk).reflec) && (isdiag(transforms(iSelChk).rot) && all(diag(transforms(iSelChk).rot) == 1))
                    transforms(iSelChk).reflec= false;
                end
                symGroup= symGroupTemplate(symGroupTemplate ~= iSelChk);
                
                cenB= mean(ptsRef);
                for sy= 1 : length(symGroup)
                    ptsSym= samplingsArray(shapeIxR).samples(segmentsArray(templateIndices(symGroup(sy))).ixSamples, :);
                    cenSym= mean(ptsSym);

                    %Find the symmetry (reflection) plane for the symmetric segment
                    tempV= (cenB-cenSym)';
                    tempV= tempV / norm(tempV);
                    maxDp= 0;
                    selAx= 0;
                    for a= 1 : 3
                        axR= zeros(3, 1);
                        axR(a)= 1;
                        dp= abs(sum(tempV .* axR));
                        if (dp > maxDp && dp > 0.95)
                           maxDp= dp;
                           selAx= a; 
                        end
                    end
                    
                    if (selAx >0)
                        transforms(symGroup(sy)).reflec= ~transforms(iSelChk).reflec;
                    else
                        transforms(symGroup(sy)).reflec= transforms(iSelChk).reflec;
                    end
                end
            end
        end
        %%%
        
        if(~doLocAlign)
            if (~err)
                err = fSavePly_FromSegments(newSegments, pathToSave, colorPerSeg);
            end
            return;
        end
        
        if (~err)
            %II. Try to optimize the novel shape, by finding contact (slot)
            %points between the neighbor segments and aligning them...
            maxScFact= 1.65;
            minScFact= 1.0;
            volThres= 0.035;
            epsilAreaF= 0.015;
            
            szShape= max(samplingsArray(shapeIxR).samples) - min(samplingsArray(shapeIxR).samples);
            volShape = szShape(1)*szShape(2)*szShape(3);
            %We use a threshold based on a percentage of the Shape's volume, as in Paper (b): 
            %this is similar to what we did for finding the neighbors of the segments
            thrDistSeg= (volShape^(1/3)) * volThres;

            contactsInfo= struct([]);
            aligned= false(size(templateIndices, 2), 1);
            adjChecked= false(size(segmentsAdj(shapeIxR).adjacencyMat));
            
            %%%
            %Organize the nodes of the Shape's graph in a list according to
            %a Breadth-First Search (References (c), (d))
            
            %First get as the initial (source) node for the search, the
            %segment having most number of neighbbors
            maxNN= 1;
            numNeigs= zeros(size(templateIndices, 2), 1);
            for i= 1 : size(templateIndices, 2)
                numN= length(find(segmentsAdj(shapeIxR).adjacencyMat(i, :) == 1));
                numNeigs(i)= numN;
                if (numN > maxNN)
                    maxNN= numN;
                    q= i;
                end
            end
                        
            sS= q;
%             sS= 1;
%             q= 1;
            while (~isempty(q))
                cIx= q(1);
                q(1)= [];
                
                ixNeighs= find(segmentsAdj(shapeIxR).adjacencyMat(cIx, :) == 1);
                
                notInS= ~ismember(ixNeighs, sS);
                sS= [sS ixNeighs(notInS)];
                q= [q ixNeighs(notInS)];
            end
            %%%
            
            for i= 1 : length(sS)
                %Find each neighbour of the new segment, according to the
                %adjacency computed for the Template shape
                ixNeighs= find(segmentsAdj(shapeIxR).adjacencyMat(sS(i), :) == 1);
                
                for nh= 1 : length(ixNeighs)
                    doCont= false;
                    
                    %Check if the segments have not been aligned
                    if (~adjChecked(sS(i), ixNeighs(nh))) || (~aligned(sS(i)) || ~aligned(ixNeighs(nh)))
                        adjChecked(sS(i), ixNeighs(nh))= true;
                        adjChecked(ixNeighs(nh), sS(i))= true;
                        
                        %%%%%
                        %For the case when we preserve some segments of the
                        %original template, no alignment will be performed
                        if (templateIndices(sS(i)) == selSegs(sS(i)))
                            aligned(sS(i))= true;
                        end
                        if ((templateIndices(ixNeighs(nh)) == selSegs(ixNeighs(nh))))
                            aligned(ixNeighs(nh))= true;
                        end
                        if (templateIndices(sS(i)) == selSegs(sS(i))) && (templateIndices(ixNeighs(nh)) == selSegs(ixNeighs(nh)))                            
                            continue;
                        end
                        %%%%%
                        
                        ptsI= newSegments(sS(i)).points;
                        ptsNH= newSegments(ixNeighs(nh)).points;

                        %Get the contact points between each pair of segments...
                        obbI= part_obb(ptsI);
                        obbNH= part_obb(ptsNH);

                        facets= obbI.faces;
                        facetsNbh= obbNH.faces;
%                         dMin= 999999;
                        dMinI= 999999;
                        dMinNH= 999999;
                        
                        %Find for each of the neighbor segments, the closest faces on the corresponding OBBs ...
                        for fsi= 1 : 6
                            %Project the points of the segment towards the plane 
                            %given by the face of the OBB
                            nI= reshape(cross(facets(fsi, 2, :)-facets(fsi, 1, :), facets(fsi, 3, :)-facets(fsi, 1, :)), [1, 3]);
                            nI= nI / norm(nI);
                            pToV1= bsxfun(@minus, ptsI, reshape(facets(fsi, 1, :), [1 3]));
                            dp= sum(bsxfun(@times, pToV1, nI), 2);
                            pProj= bsxfun(@minus, ptsI, bsxfun(@times, dp, nI));

                            %Get the sum of the distances between the projected points 
                            %and the corresponding closest point on the neighbor segment
                            [~, distBS]= knnsearch(ptsNH, pProj);
                            tempD= sum(distBS);
                            if (tempD < dMinI)
                                dMinI= tempD;
                                minFI= fsi;
                            end
                        end
                                                                                        
                        for fsn= 1 : 6
                            %Project the points of the segment towards the plane 
                            %given by the face of the OBB
                            nNh= reshape(cross(facetsNbh(fsn, 2, :)-facetsNbh(fsn, 1, :), facetsNbh(fsn, 3, :)-facetsNbh(fsn, 1, :)), [1, 3]);
                            nNh= nNh / norm(nNh);
                            pToV1= bsxfun(@minus, ptsNH, reshape(facetsNbh(fsn, 1, :), [1 3]));
                            dp= sum(bsxfun(@times, pToV1, nNh), 2);
                            pProj= bsxfun(@minus, ptsNH, bsxfun(@times, dp, nNh));

                            %Get the sum of the distances between the projected points 
                            %and the corresponding closest point on the neighbor segment
                            [~, distBS]= knnsearch(ptsI, pProj);
                            tempD= sum(distBS);
                            if (tempD < dMinNH)
                               dMinNH= tempD;
                               minFNh= fsn;
                            end
                        end
                            
                        %Get the four corners of the closest faces of each OBB of the new segments
                        cornFI= reshape(facets(minFI, :, :), [4 3]);
                        cornFNH= reshape(facetsNbh(minFNh, :, :), [4 3]);
                        
                        %Get area of the closest faces on each OBB
                        areaCfI= norm(cornFI(1, :) - cornFI(2, :)) * norm(cornFI(3, :) - cornFI(4, :));
                        areaCfNh= norm(cornFNH(1, :) - cornFNH(2, :)) * norm(cornFNH(3, :) - cornFNH(4, :));
                        
                        %If the area of the closest faces is similar, we
                        %use the point closest to each corner as contact
                        ptsClose= zeros(4, 1);
                        if (abs(areaCfI - areaCfNh) < epsilAreaF)
                            closestI3= cornFI;
                            closestNH3= cornFNH;
                        else
                            %If area of each face is not similar, check which face has the smallest area
                            if (areaCfI < areaCfNh)
                                closestI3= cornFI;
                                
                                %Project the corners of segment "i" towards the closest 
                                %face of the OBB of segment "nh"
                                nNh= reshape(cross(facetsNbh(minFNh, 2, :)-facetsNbh(minFNh, 1, :), facetsNbh(minFNh, 3, :)-facetsNbh(minFNh, 1, :)), [1, 3]);
                                nNh= nNh / norm(nNh);
                                pToV1= bsxfun(@minus, cornFI, reshape(facetsNbh(minFNh, 1, :), [1 3]));
                                dp= sum(bsxfun(@times, pToV1, nNh), 2);
                                pPrj= bsxfun(@minus, cornFI, bsxfun(@times, dp, nNh));
                                
                                tmpIxClose= knnsearch(ptsNH, pPrj, 'K', 4);
                                ptsClose(1)= tmpIxClose(1);
                                for cluni= 2 : 4
                                    for theK= 1 : 4
                                        %We get unique points: we want 4 distinct contacts
                                        if (~ismember(tmpIxClose(cluni, theK), ptsClose))
                                            ptsClose(cluni)= tmpIxClose(cluni, theK);
                                            break;
                                        end
                                    end
                                end
                                closestNH3= ptsNH(ptsClose, :);
                            else
                                closestNH3= cornFNH;
                                
                                %Project the corners of segment "nh" towards the closest 
                                %face of the OBB of segment "i"
                                nI= reshape(cross(facets(minFI, 2, :)-facets(minFI, 1, :), facets(minFI, 3, :)-facets(minFI, 1, :)), [1, 3]);
                                nI= nI / norm(nI);
                                pToV1= bsxfun(@minus, cornFNH, reshape(facets(minFI, 1, :), [1 3]));
                                dp= sum(bsxfun(@times, pToV1, nI), 2);
                                pPrj= bsxfun(@minus, cornFNH, bsxfun(@times, dp, nI));
                                
                                tmpIxClose= knnsearch(ptsI, pPrj, 'K', 4);
                                ptsClose(1)= tmpIxClose(1);
                                for cluni= 2 : 4
                                    for theK= 1 : 4
                                        %We get unique points: we want 4 distinct contacts
                                        if (~ismember(tmpIxClose(cluni, theK), ptsClose))
                                            ptsClose(cluni)= tmpIxClose(cluni, theK);
                                            break;
                                        end
                                    end
                                end
                                closestI3= ptsI(ptsClose, :);
                            end
                        end

                        closestI= closestI3;
                        closestNH= closestNH3;

                        %As other option to see if the alignment is
                        %required, we check the average distance 
                        %between the candidate contact points
                        allDs= pdist2(closestI, closestNH);
                        dbtwS= diag(allDs);
                        doCont= (mean(dbtwS) > thrDistSeg);
                        sta= 0;
                        if (doCont)
                            contactsInfo(sS(i), ixNeighs(nh)).points= closestI;
                            contactsInfo(ixNeighs(nh), sS(i)).points= closestNH;
                            if (templateIndices(sS(i)) == selSegs(sS(i)))
                                %If one of the segments is from the template, we will apply the transformation to the other segment
                                %NOTE: if both adjacent segments are from the template, this code will not be
                                %executed: see the conditions (if) above
                                srcPoints= contactsInfo(ixNeighs(nh), sS(i)).points;
                                targPoints= contactsInfo(sS(i), ixNeighs(nh)).points;
                                sta= ixNeighs(nh);
                            elseif (templateIndices(ixNeighs(nh)) == selSegs(ixNeighs(nh)))
                                %If one of the segments is from the template, we will apply the transformation to the other segment
                                %NOTE: if both adjacent segments are from the template, this code will not be
                                %executed: see the conditions (if) above
                                srcPoints= contactsInfo(sS(i), ixNeighs(nh)).points;
                                targPoints= contactsInfo(ixNeighs(nh), sS(i)).points;
                                sta= sS(i);
                            else
                                if (~aligned(sS(i)) && ~aligned(ixNeighs(nh)))    
                                    szI= max(ptsI) - min(ptsI);
                                    szNH= max(ptsNH) - min(ptsNH);
                                    %At first, none of the two segments has been aligned, 
                                    %so we select as target the segment with greater volume
                                    volI= szI(1)*szI(2)*szI(3);
                                    volNH= szNH(1)*szNH(2)*szNH(3);
                                    if (volI < volNH)
                                        srcPoints= contactsInfo(sS(i), ixNeighs(nh)).points;
                                        targPoints= contactsInfo(ixNeighs(nh), sS(i)).points;
                                        sta= sS(i);
                                    else    
                                        srcPoints= contactsInfo(ixNeighs(nh), sS(i)).points;
                                        targPoints= contactsInfo(sS(i), ixNeighs(nh)).points;
                                        sta= ixNeighs(nh);
                                    end
                                else
                                    %If one of the two segments has been aligned, we'll select 
                                    %as target the segment that was already aligned
                                    if (aligned(sS(i)))
                                        srcPoints= contactsInfo(ixNeighs(nh), sS(i)).points;
                                        targPoints= contactsInfo(sS(i), ixNeighs(nh)).points;
                                        sta= ixNeighs(nh);
                                    elseif (aligned(ixNeighs(nh)))
                                        srcPoints= contactsInfo(sS(i), ixNeighs(nh)).points;
                                        targPoints= contactsInfo(ixNeighs(nh), sS(i)).points;
                                        sta= sS(i);
                                    end
                                end 
                            end
                            
                            if (sta > 0)
                                mainAx= segmentsArray(templateIndices(sta)).segmentMainAx;
                            
                                if ((sta == sS(i)) && (mainAx ~= segmentsArray(templateIndices(ixNeighs(nh))).segmentMainAx)) || ...
                                        ((sta == ixNeighs(nh)) && (mainAx ~= segmentsArray(templateIndices(sS(i))).segmentMainAx))
                                    if (all(abs(closestI(:, mainAx) - closestNH(:, mainAx)) < thrDistSeg))
                                        sta= 0;
                                    end
                                end
                            end
                        end
                      
                        if (sta > 0)
                            %Next, we adapt the code for alignment (translation and scaling) between two 
                            %sets of points by Han Liu, which in turn, is based on code from Paper (a)
                            numCl= size(srcPoints, 1);
                            centSrc= mean(srcPoints);
                            centTarg= mean(targPoints);

                            newScX= 0;
                            newScY= 0;
                            newScZ= 0;
                            for ixfs= 1 : numCl
                                pS= srcPoints(ixfs) - centSrc;
                                pT= targPoints(ixfs) - centTarg;
                                newScX = newScX + (pT(1) / pS(1));
                                newScY = newScY + (pT(2) / pS(2));
                                newScZ = newScZ + (pT(3) / pS(3));
                            end
                            newScX = newScX / numCl;
                            newScY = newScY / numCl;
                            newScZ = newScZ / numCl;

                            newScX(isnan(newScX))= minScFact;
                            newScY(isnan(newScY))= minScFact;
                            newScZ(isnan(newScZ))= minScFact;

                            %%%
                            %We try another scaling factor, and we select the largest
                            %valid factor. This is because, in general, we're trying to 
                            %fill any hole between neighbor segments in our shape
                            %NOTE: this additional factor and all the additional heuristics below, 
                            %are Not from Han's original code
                            sz2= (max(targPoints)-min(targPoints)) ./ (max(srcPoints)-min(srcPoints));

                            if (sz2(1)>newScX) && (sz2(1)<maxScFact)
                                newScX= sz2(1);
                            end
                            if (sz2(2)>newScY) && (sz2(2)<maxScFact)
                                newScY= sz2(2);
                            end
                            if (sz2(3)>newScZ) && (sz2(3)<maxScFact)
                                newScZ= sz2(3);
                            end

                            newScX(newScX < minScFact)= minScFact;
                            newScY(newScY < minScFact)= minScFact;
                            newScZ(newScZ < minScFact)= minScFact;

                            newScX(newScX > maxScFact)= maxScFact;
                            newScY(newScY > maxScFact)= maxScFact;
                            newScZ(newScZ > maxScFact)= maxScFact;

                            newmScale= [newScX 0 0; 0 newScY 0; 0 0 newScZ];
                            
                            ptsOfTemplate= samplingsArray(shapeIxR).samples(segmentsArray(templateIndices(sta)).ixSamples, :);
                            szTemplate= max(ptsOfTemplate)-min(ptsOfTemplate);

                            %Try the new scaling, and measure vs. the original segment of the Template, 
                            %to see how well the new segment fits
                            ptsCurrent= newSegments(sta).points;
                            ptsCurrent= ptsCurrent * newmScale;
                            szNew= max(ptsCurrent)-min(ptsCurrent);

                            if (szNew(1) > szTemplate(1))
                                newScX= newScX * (1-(szNew(1)-szTemplate(1)));
                                newScX(newScX < 1.0)= 1.0;
                            end
                            if (szNew(2) > szTemplate(2))
                                newScY= newScY * (1-(szNew(2)-szTemplate(2)));
                                newScY(newScY < 1.0)= 1.0;
                            end
                            if (szNew(3) > szTemplate(3))
                                newScZ= newScZ * (1-(szNew(3)-szTemplate(3)));
                                newScZ(newScZ < 1.0)= 1.0;
                            end

                            if (neighIsSymm(sta))
                                %For the special case of segments that are neighbours and also symmetric, 
                                %we do not require any scaling through the Y coord (upright vector)
                                newScY= 1;
                            end
                            newmScale= [newScX 0 0; 0 newScY 0; 0 0 newScZ];
                            %%%
                            
                            %Save current centroid of segment
                            currMean= mean(newSegments(sta).points);
                            
                            [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(sta), segmentsArray, templateIndices);
                            if (err)
                                break;
                            end
                            
                            newSegments(sta).points= newSegments(sta).points * newmScale;

                            %Recenter the segment after the scaling
                            deltaNLoc= currMean - mean(newSegments(sta).points);
                            newSegments(sta).points= newSegments(sta).points + repmat(deltaNLoc, size(newSegments(sta).points, 1), 1);
                            %Save the new scaling
                            transforms(sta).scale= transforms(sta).scale * newmScale;
                            
                            aligned(sta)= true;
                            
                            %%%
                            %We apply a similar transformation for the
                            %symmetric segment...
                            
                            %Calculate and save the translation that was applied
                            newSegTr= segmentsArray(selSegs(sta));
                            ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                            ptsNewTr= ptsNewTr * transforms(sta).rot;
                            ptsNewTr= ptsNewTr * transforms(sta).scale;
                            deltaP= mean(newSegments(sta).points) - mean(ptsNewTr);
                            transforms(sta).translat= deltaP;
                            
                            if (~isempty(symGroupTemplate))
                                symGroup= symGroupTemplate(symGroupTemplate ~= sta);
                                
                                
                                for sy= 1 : length(symGroup)
                                    %Save current centroid of segment
                                    currMean= mean(newSegments(symGroup(sy)).points);                                    
                                    newSegments(symGroup(sy)).points= newSegments(symGroup(sy)).points * newmScale;

                                    %Recenter the segment after the scaling
                                    deltaNLoc= currMean - mean(newSegments(symGroup(sy)).points);
                                    newSegments(symGroup(sy)).points= newSegments(symGroup(sy)).points + ...
                                        repmat(deltaNLoc, size(newSegments(symGroup(sy)).points, 1), 1);
                                    transforms(symGroup(sy)).scale= transforms(sta).scale;
                                    
                                    %Calculate and save the translation that was applied
                                    ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                                    ptsNewTr= ptsNewTr * transforms(symGroup(sy)).rot;
                                    ptsNewTr= ptsNewTr * transforms(symGroup(sy)).scale;

                                    deltaP= mean(newSegments(symGroup(sy)).points) - mean(ptsNewTr);
                                    transforms(symGroup(sy)).translat= deltaP;
                                end
                                aligned(symGroup)= true;                                
                            end
                            %%%
                        elseif (templateIndices(sS(i)) ~= selSegs(sS(i)))
                            %For cases where we did not apply the additional alignment...
                            [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(sS(i)), segmentsArray, templateIndices);
                            if (err)
                                break;
                            end
                            
                            %Calculate and save the translation that was applied
                            newSegTr= segmentsArray(selSegs(sS(i)));
                            ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                            ptsNewTr= ptsNewTr * transforms(sS(i)).rot;
                            ptsNewTr= ptsNewTr * transforms(sS(i)).scale;
                            deltaP= mean(newSegments(sS(i)).points) - mean(ptsNewTr);
                            transforms(sS(i)).translat= deltaP;
                            
                            if (~isempty(symGroupTemplate))
                                symGroup= symGroupTemplate(symGroupTemplate ~= sS(i));
                                
                                for sy= 1 : length(symGroup)
                                    %Calculate and save the translation that was applied
                                    ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                                    ptsNewTr= ptsNewTr * transforms(symGroup(sy)).rot;
                                    ptsNewTr= ptsNewTr * transforms(symGroup(sy)).scale;

                                    deltaP= mean(newSegments(symGroup(sy)).points) - mean(ptsNewTr);
                                    transforms(symGroup(sy)).translat= deltaP;
                                end
                                
                                if (aligned(sS(i)))
                                    aligned(symGroup)= true;
                                end
                            end
                        end
                    else
                        if (~any(transforms(sS(i)).translat)) && (templateIndices(sS(i)) ~= selSegs(sS(i)))
                            %For cases where we did not apply the additional alignment...
                            [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(sS(i)), segmentsArray, templateIndices);
                            if (err)
                                break;
                            end
                            
                            %Calculate and save the translation that was applied
                            newSegTr= segmentsArray(selSegs(sS(i)));
                            ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                            ptsNewTr= ptsNewTr * transforms(sS(i)).rot;
                            ptsNewTr= ptsNewTr * transforms(sS(i)).scale;
                            deltaP= mean(newSegments(sS(i)).points) - mean(ptsNewTr);
                            transforms(sS(i)).translat= deltaP;
                            
                            if (~isempty(symGroupTemplate))
                                symGroup= symGroupTemplate(symGroupTemplate ~= sS(i));
                                
                                for sy= 1 : length(symGroup)
                                    %Calculate and save the translation that was applied
                                    ptsNewTr= samplingsArray(newSegTr.ixShape).samples(newSegTr.ixSamples, :);
                                    ptsNewTr= ptsNewTr * transforms(symGroup(sy)).rot;
                                    ptsNewTr= ptsNewTr * transforms(symGroup(sy)).scale;

                                    deltaP= mean(newSegments(symGroup(sy)).points) - mean(ptsNewTr);
                                    transforms(symGroup(sy)).translat= deltaP;
                                end
                                
                                if (aligned(sS(i)))
                                    aligned(symGroup)= true;
                                end
                            end
                        end
                    end
                end
                
                if (err)
                    break;
                end
            end            
        end

        if (~err)
            err = fSavePly_FromSegments(newSegments, pathToSave, colorPerSeg);
        end
    catch ME
        err= -1;
        sErr= ['Error synthesizing a new shape, (error in fSynthesizeNovelShape): ' ME.message];
        errordlg(sErr);
    end
end

function [ptsNew, mScale]= fTransAndScalePoints(ptsIni, ptsRef, szRef)
    %We scale the new segment according to the proportions of the template's segment
    szNew= max(ptsIni)-min(ptsIni);

    scaleFactX= szRef(1) / szNew(1);
    scaleFactY= szRef(2) / szNew(2);
    scaleFactZ= szRef(3) / szNew(3);

    mScale= [scaleFactX 0 0; 0 scaleFactY 0; 0 0 scaleFactZ];

    ptsNew= ptsIni * mScale;

    %We translate the new segment towards the position of the template's segment
    deltaP= mean(ptsRef) - mean(ptsNew);
%                 deltaP= min(ptsRef) - min(ptsNew);
    ptsNew= ptsNew + repmat(deltaP, size(ptsNew, 1), 1);

    delta2= min(ptsRef) - min(ptsNew);
    ptsNew= ptsNew + repmat(delta2, size(ptsNew, 1), 1);
end