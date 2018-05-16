%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Aug/2017 - XXX/2017
%Funcion: fSynthesizeMeshFromPoints
%Input:
%   meshesArray:        array of structures containing the Meshes
%   segmentsArray:      structure array containing the segments that were previously computed for all 
%                       the input shapes
%   samplingsArray:     structure array containing the samples of the original Shapes
%   selSegs:            indices of the segments that were selected as replacements for each segment 
%                       of the template, (see fSampleSegsForNewShape)
%   pathToSave:         full path (including file name and extension) for
%                       saving the new mesh, in an OBJ file
%Output:
%   err:        -1 if some error ocurrs; 0, otherwise
%   newMesh:    structure containing the infor of the new segments used for the synthesized shape
%
% Function that synthesizes a new shape using the new segments ("selSegs") to replace 
% each segment of a Template Shape, and saves the result to a PLY file
%
% (See also mesh_get_con, fSynthesizeNovelShape)
%%%%%
%(1) References: 
%   (a) Zheng, Y., et al. "Smart Variations: Functional Substructures for Part Compatibility" (2013)
%   (b) Prakhar, J., et al. "Assembly-based conceptual 3D modeling with unlabeled components using probabilistic factor graph" (2016)
%   (c) https://en.wikipedia.org/wiki/Breadth-first_search
%%%%%
%%%
%--------------------------------------------------------------------------
% function [err, newMesh] = fSynthesizeMeshFromPoints(meshSegments, segmentsArray, selSegs, transforms, templateIndices, segmentsAdj, doLocAlign)
% function [err, newMesh] = fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSegs, transforms, templateIndices, segmentsAdj, doLocAlign)
function [err, newMesh] = fSynthesizeMeshFromPointsNS(meshSegments, segmentsArray, samplingsArray, selSegs, transforms, templateIndices, segmentsAdj, ...
    doLocAlign, savePath, doSaveNoAlign)
    err= 0;
    try
        if length(selSegs) ~= length(templateIndices)
            sErr= 'Error: the number of segments of the template must match the number of the new segments. Check the input parameters';
            errordlg(sErr);
            return;
        end
        if (~exist('doLocAlign', 'var'))
            doLocAlign= true;
        end
        if (~exist('doSaveNoAlign', 'var'))
            doSaveNoAlign= false;
        end
        
        shapeIxR= segmentsArray(templateIndices(1)).ixShape;
%         volThres= 0.01;
% %         szShape= max(meshesArray(shapeIxR).vertices) - min(meshesArray(shapeIxR).vertices);
%         szShape= max(samplingsArray(shapeIxR).samples) - min(samplingsArray(shapeIxR).samples);
%         volShape = szShape(1)*szShape(2)*szShape(3);
% %         %We use a threshold based on a percentage of the Shape's volume, as in Paper (b): 
% %         %this is similar to what we did for finding the neighbors of the segments
%         thrDistSeg= (volShape^(1/3)) * volThres;
%         thrDistSeg= 0.0025;
        epsil1= 0.0025;
        
        newMesh= [];
        vers= [];
        tris= [];
        ixNextIx= 0;
        meshParts= [];
        %STEP 1: go back to the meshes and extract the faces that belong to each segment of the new shape. 
        %To create the new mesh, we transform each original mesh part
        %according to the transformation that we computed before for the points        
        for i= 1 : length(selSegs)
            meshOfSeg= meshSegments(selSegs(i));
            
            %Apply the transformation to the mesh vertices
            vxs= meshOfSeg.vertices;
            vxs= vxs * transforms(i).rot;
            vxs= vxs * transforms(i).scale;

            vxs= vxs + repmat(transforms(i).translat, size(vxs, 1), 1);
                        
            if (~isempty(vers))
                iniLenVers= length(vers) + 1;
                vers= cat(1, vers, vxs);
            else
                iniLenVers= 1;
                vers= vxs;
            end
            
            %Get the face indices for the new mesh
%             [~, ixF]= ismember(meshOfSeg.faces(facesPerSeg, :), vIx);
            ixF= meshOfSeg.faces;
            
            %If the transformation implied a reflection, change the order of the vertices, 
            %to avoid having faces with wrong orientation
            if (transforms(i).reflec)
                %We exchange vertices 1 and 3, (we could exchange any pair)
                tmpF= ixF(:, 3);
                ixF(:, 3)= ixF(:, 1);
                ixF(:, 1)= tmpF;
            end
            
            ixF= ixF + ixNextIx;
            if (~isempty(tris))
                tris= cat(1, tris, ixF);
            else
                tris= ixF;
            end
            ixNextIx= ixNextIx + size(vxs, 1);
            
            mP.vertices= vxs;
            mP.faces= ixF;
            mP.ixVecIni= iniLenVers;
            mP.ixVecEnd= length(vers);
            meshParts= [meshParts mP];
        end
        newMesh.vertices = vers;
        newMesh.faces= tris;
        
        if (exist('savePath', 'var'))
            if (~isempty(savePath) && doSaveNoAlign)
                [fPath, fName, fExt]= fileparts(savePath);                
                mesh_write_obj(newMesh, fullfile(fPath, [fName '-NoA.' fExt]));
            end
        end
        
        if (doLocAlign)            
            neighIsSymm= false(length(templateIndices), 1);
            
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
                
                %Also, save a variable saying if a neighbor is also symmetric
                ixNeighs= find(segmentsAdj(shapeIxR).adjacencyMat(i, :) == 1);
                [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(i), segmentsArray, templateIndices);
                refSymIsNeigh= find(ismember(symGroupTemplate, ixNeighs), 1);
                if (~isempty(refSymIsNeigh))
                    neighIsSymm(i)= true;
                end
            end
            
            %%%
            %Organize the nodes of the Shape's graph in a list according to
            %a Breadth-First Search (Reference (c))
            sS= q;
            while (~isempty(q))
                cIx= q(1);
                q(1)= [];

                ixNeighs= find(segmentsAdj(shapeIxR).adjacencyMat(cIx, :) == 1);                
                [~, ixTmp]= sort(numNeigs(ixNeighs), 'descend');
                ixNeighs= ixNeighs(ixTmp);
                
                notInS= ~ismember(ixNeighs, sS);
                sS= [sS ixNeighs(notInS)];
                q= [q ixNeighs(notInS)];
            end
            %%%

            aligned= false(size(templateIndices, 2), 1);
            adjChecked= false(size(segmentsAdj(shapeIxR).adjacencyMat));
            contactsInfo= struct([]);

            for i= 1 : length(sS)
                %Find each neighbour of the new segment, according to the
                %adjacency computed for the Template shape
                ixNeighs= find(segmentsAdj(shapeIxR).adjacencyMat(sS(i), :) == 1);
                for nh= 1 : length(ixNeighs)                    

                    ixIniP= meshParts(sS(i)).ixVecIni;
                    ixEndP= meshParts(sS(i)).ixVecEnd;
                    ptsI= newMesh.vertices(ixIniP:ixEndP, :);

                    ixIniP= meshParts(ixNeighs(nh)).ixVecIni;
                    ixEndP= meshParts(ixNeighs(nh)).ixVecEnd;
                    ptsNH= newMesh.vertices(ixIniP:ixEndP, :);

                    doCont= false;
                    %%%
                    if (segmentsAdj(shapeIxR).doAlign(sS(i), ixNeighs(nh)) == 0)
                        adjChecked(sS(i), ixNeighs(nh))= true;
                        adjChecked(ixNeighs(nh), sS(i))= true;
                    end
                    %%%

                    if ((~adjChecked(sS(i), ixNeighs(nh))) || (~adjChecked(ixNeighs(nh), sS(i)))) && (~aligned(sS(i)) || ~aligned(ixNeighs(nh)))
                        contactsInfo(sS(i), ixNeighs(nh)).points= [];
                        contactsInfo(ixNeighs(nh), sS(i)).points= [];

                        contactsInfo(sS(i), ixNeighs(nh)).ixSamp= [];
                        contactsInfo(ixNeighs(nh), sS(i)).ixSamp= [];

                        bbI= fGetBoundBoxForPoints(ptsI);
                        bbNH= fGetBoundBoxForPoints(ptsNH);
                        facetsI= bbI.faces;
                        facetsNh= bbNH.faces;

                        szI= max(ptsI) - min(ptsI);
                        szNH= max(ptsNH) - min(ptsNH);
                        volI= szI(1)*szI(2)*szI(3);
                        volNH= szNH(1)*szNH(2)*szNH(3);

                        %Same main axis of orientation?
                        if (segmentsArray(templateIndices(sS(i))).segmentMainAx == segmentsArray(templateIndices(ixNeighs(nh))).segmentMainAx)
                            mnAx= segmentsArray(templateIndices(sS(i))).segmentMainAx;
                            [clFI, clFNh]= fGetClosesOBBFaces(mnAx, facetsI, facetsNh);

                            %Get corners of each closest face
                            cornFI= reshape(facetsI(clFI, :, :), [4 3]);
                            cornFNH= reshape(facetsNh(clFNh, :, :), [4 3]);

                            ptsClose= zeros(4, 1);
                            tmpIxClose= knnsearch(ptsI, cornFI, 'K', 4);
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
                            contactsInfo(sS(i), ixNeighs(nh)).points= ptsI(ptsClose, :);
                            contactsInfo(sS(i), ixNeighs(nh)).ixSamp= ptsClose;

                            tmpIxClose= knnsearch(ptsNH, cornFNH, 'K', 4);
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
                            contactsInfo(ixNeighs(nh), sS(i)).points= ptsNH(ptsClose, :);
                            contactsInfo(ixNeighs(nh), sS(i)).ixSamp= ptsClose;

                            doCont= all(abs(contactsInfo(sS(i), ixNeighs(nh)).points(:,mnAx) - contactsInfo(ixNeighs(nh), sS(i)).points(:,mnAx)) < 0.25);
%                         end
                        else
                            if (neighIsSymm(sS(i))) || (neighIsSymm(ixNeighs(nh)))
                                %If we have one of the segments as "neighbour and symmetric", always use the other
                                %segment's axis as main axis
                                if (neighIsSymm(sS(i)))
                                    mnAx= segmentsArray(templateIndices(ixNeighs(nh))).segmentMainAx;
                                elseif (neighIsSymm(ixNeighs(nh)))
                                    mnAx= segmentsArray(templateIndices(sS(i))).segmentMainAx;
                                end

                                [clFI, clFNh]= fGetClosesOBBFaces(mnAx, facetsI, facetsNh);

                                cornFI= reshape(facetsI(clFI, :, :), [4 3]);
                                cornFNH= reshape(facetsNh(clFNh, :, :), [4 3]);
                            else
                                %Check the main axis of each segment and keep as closest faces the ones that have
                                %the min. total distance between the corners of the face and the points of the
                                %other segment
                                [clFI, clFNh]= fGetClosesOBBFaces(segmentsArray(templateIndices(sS(i))).segmentMainAx, facetsI, facetsNh);
                                cornFI1= reshape(facetsI(clFI, :, :), [4 3]);
                                cornFNH1= reshape(facetsNh(clFNh, :, :), [4 3]);
                                
                                [~, dMinI1]= knnsearch(ptsI, cornFNH1);
                                [~, dMinNh1]= knnsearch(ptsNH, cornFI1);
                                dMin1= mean(cat(1, dMinI1, dMinNh1));
                                
                                [clFI, clFNh]= fGetClosesOBBFaces(segmentsArray(templateIndices(ixNeighs(nh))).segmentMainAx, facetsI, facetsNh);
                                cornFI2= reshape(facetsI(clFI, :, :), [4 3]);
                                cornFNH2= reshape(facetsNh(clFNh, :, :), [4 3]);
                                
                                [~, dMinI2]= knnsearch(ptsI, cornFNH2);
                                [~, dMinNh2]= knnsearch(ptsNH, cornFI2);
                                dMin2= mean(cat(1, dMinI2, dMinNh2));
                                
                                doCont= (dMin1 > 0.03 & dMin2 > 0.03) & (dMin1 < 0.15 & dMin2 < 0.15);
                                if (doCont)
                                    if (dMin1 < dMin2)
                                        cornFI= cornFI1;
                                        cornFNH= cornFNH1;
                                    else
                                        cornFI= cornFI2;
                                        cornFNH= cornFNH2;
                                    end
                                end
                            end
                            
                            if (doCont)
                                %Get area of the closest faces on each OBB
                                areaCfI= norm(cornFI(1, :) - cornFI(2, :)) * norm(cornFI(3, :) - cornFI(4, :));
                                areaCfNh= norm(cornFNH(1, :) - cornFNH(2, :)) * norm(cornFNH(3, :) - cornFNH(4, :));

                                if (areaCfI < areaCfNh)
                                    ptsClose= zeros(4, 1);
                                    tmpIxClose= knnsearch(ptsI, cornFI, 'K', 4);
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
                                    contactsInfo(sS(i), ixNeighs(nh)).points= ptsI(ptsClose, :);
                                    contactsInfo(sS(i), ixNeighs(nh)).ixSamp= ptsClose;

                                    tmpIxClose= knnsearch(ptsNH, contactsInfo(sS(i), ixNeighs(nh)).points, 'K', 4);
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
                                    contactsInfo(ixNeighs(nh), sS(i)).points= ptsNH(ptsClose, :);
                                    contactsInfo(ixNeighs(nh), sS(i)).ixSamp= ptsClose;
                                else
                                    tmpIxClose= knnsearch(ptsNH, cornFNH, 'K', 4);
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
                                    contactsInfo(ixNeighs(nh), sS(i)).points= ptsNH(ptsClose, :);
                                    contactsInfo(ixNeighs(nh), sS(i)).ixSamp= ptsClose;

                                    tmpIxClose= knnsearch(ptsI, contactsInfo(ixNeighs(nh), sS(i)).points, 'K', 4);
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
                                    contactsInfo(sS(i), ixNeighs(nh)).points= ptsI(ptsClose, :);
                                    contactsInfo(sS(i), ixNeighs(nh)).ixSamp= ptsClose;
                                end
                            end
                        end
                        
                        if (doCont)
                            minPS= min(ptsI) + epsil1;
                            maxPS= max(ptsI) - epsil1;
                            doCont= ~(all(all(bsxfun(@ge, contactsInfo(ixNeighs(nh), sS(i)).points, minPS), 2) & ...
                                all(bsxfun(@le, contactsInfo(ixNeighs(nh), sS(i)).points, maxPS), 2)));
                        end
                        if (doCont)
                            minPS= min(ptsNH) + epsil1;
                            maxPS= max(ptsNH) - epsil1;
                            doCont= ~(all(all(bsxfun(@ge, contactsInfo(sS(i), ixNeighs(nh)).points, minPS), 2) & ...
                                all(bsxfun(@le, contactsInfo(sS(i), ixNeighs(nh)).points, maxPS), 2)));
                        end

                        sta= 0;
                        if (doCont)
                            if (~aligned(sS(i)) && ~aligned(ixNeighs(nh))) || (aligned(sS(i)) && aligned(ixNeighs(nh)))
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

                        if (doCont && sta > 0)
                            if (~neighIsSymm(sta))
                                [newmScale, transV]= fCalcLocalAlignment(srcPoints, targPoints, 1.45, 0.8);

                                if (~isempty(newmScale))
                                    %%%
                                    mainAx= segmentsArray(templateIndices(sta)).segmentMainAx;
                                    if (newmScale(mainAx, mainAx) < 1.0)
                                        newmScale(mainAx, mainAx)= 1.0;
                                    end
                                    %%%

                                    ixIniP= meshParts(sta).ixVecIni;
                                    ixEndP= meshParts(sta).ixVecEnd;

                                    currMean= mean(newMesh.vertices(ixIniP:ixEndP, :));

                                    newMesh.vertices(ixIniP:ixEndP, :)= newMesh.vertices(ixIniP:ixEndP, :) * newmScale;
                                    %After scaling, we recenter the segment
                                    newDeltaP= currMean - mean(newMesh.vertices(ixIniP:ixEndP, :));
                                    newMesh.vertices(ixIniP:ixEndP, :)= newMesh.vertices(ixIniP:ixEndP, :) + ...
                                        repmat(newDeltaP, size(newMesh.vertices(ixIniP:ixEndP, :), 1), 1);

                                    %Translate the segment
                                    %%%
                                    if (abs(transV(mainAx)) > 0.025)
                                        transV(mainAx)= sign(transV(mainAx)) * 0.025;
                                    end

                                    [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(sta), segmentsArray, templateIndices);
                                    if isempty(symGroupTemplate)
                                        transV(3)= 0.0;
                                    end
                                    newMesh.vertices(ixIniP:ixEndP, :)= newMesh.vertices(ixIniP:ixEndP, :) + ...
                                            repmat(transV, size(newMesh.vertices(ixIniP:ixEndP, :), 1), 1);

                                    %%%
                                    aligned(sta)= true;
                                end
                            end
                        end
                    end

                    adjChecked(sS(i), ixNeighs(nh))= true;
                    adjChecked(ixNeighs(nh), sS(i))= true;
                end
            end
            
            %%%
            %After the alignment, check the symmetry
            segDone= false(length(templateIndices), 1);
            for i= 1 : length(templateIndices)
                if (~segDone(i))
                    %Find the label of the symmetric segments w.r.t the current template's segment
                    [err, symGroupTemplate] = fGetSymmetryChain(templateIndices(i), segmentsArray, templateIndices);
                    if (err)
                        break;
                    end

                    if (~isempty(symGroupTemplate))
                        meshOrig= meshSegments(templateIndices(i));
                        
                        ixIniP= meshParts(i).ixVecIni;
                        ixEndP= meshParts(i).ixVecEnd;
                        ptsI= newMesh.vertices(ixIniP:ixEndP, :);

                        %Get size of current template segment
                        vxs= meshOrig.vertices;
                        szOrigI= max(vxs) - min(vxs);
                                                
                        %Get size of current new segment
                        szI= max(ptsI) - min(ptsI);

                        %As a simple heuristic, take as "better" segment, the one 
                        %whose size is more similar to the original template's segment
                        szTarg= [];
                        symGroup= symGroupTemplate(symGroupTemplate ~= i);
                        for sy= 1 : length(symGroup)
                            ixIniP= meshParts(symGroup(sy)).ixVecIni;
                            ixEndP= meshParts(symGroup(sy)).ixVecEnd;
                            ptsSym= newMesh.vertices(ixIniP:ixEndP, :);
                            szSym= max(ptsSym) - min(ptsSym);

                            if any(szI ~= szSym)
                                if (abs(szOrigI - szI) < abs(szOrigI - szSym))
                                    szTarg= szI;                                    
                                else
                                    szTarg= szSym;
                                end
                            end                            
                        end
                        
                        if (~isempty(szTarg))
                            for sy= 1 : length(symGroupTemplate)
                                ixIniP= meshParts(symGroupTemplate(sy)).ixVecIni;
                                ixEndP= meshParts(symGroupTemplate(sy)).ixVecEnd;
                                ptsSym= newMesh.vertices(ixIniP:ixEndP, :);
                                szSym= max(ptsSym) - min(ptsSym);

                                scaleFactX= szTarg(1) / szSym(1);
                                scaleFactY= szTarg(2) / szSym(2);
                                scaleFactZ= szTarg(3) / szSym(3);

                                mScale= [scaleFactX 0 0; 0 scaleFactY 0; 0 0 scaleFactZ];

                                currMean= mean(newMesh.vertices(ixIniP:ixEndP, :));

                                newMesh.vertices(ixIniP:ixEndP, :)= newMesh.vertices(ixIniP:ixEndP, :) * mScale;

                                %After scaling, we recenter the segment
                                newDeltaP= currMean - mean(newMesh.vertices(ixIniP:ixEndP, :));
                                newMesh.vertices(ixIniP:ixEndP, :)= newMesh.vertices(ixIniP:ixEndP, :) + ...
                                    repmat(newDeltaP, size(newMesh.vertices(ixIniP:ixEndP, :), 1), 1);
                            end
                        end
                        segDone(symGroupTemplate)= true;
                    end
                end
            end
            %%%
        
            if (exist('savePath', 'var'))
                if (~isempty(savePath))
                    mesh_write_obj(newMesh, savePath);
                end
            end
        end
        
    catch ME
        sErr= ['Error synthesizing a new shape, (error in fSynthesizeMeshFromPoints): ' ME.message];
        errordlg(sErr);
        err= -1;
    end
end

% function [face1, face2] = fGetClosesOBBFaces(segMainAx, facets1, facets2)
%     dMin= 999999;
%     if (segMainAx == 1)
%         %In this case, by construction of the Bound. Box, we know that the closest faces are 1 and 2 or viceversa
%         fI1= reshape(facets1(1, :, :), [4, 3]);
%         fI2= reshape(facets1(2, :, :), [4, 3]);
% 
%         fNh1= reshape(facets2(1, :, :), [4, 3]);
%         fNh2= reshape(facets2(2, :, :), [4, 3]);
% 
%         dBts= diag(pdist2(fI1, fNh2));
%         if (mean(dBts) < dMin)
%             face1= 1;
%             face2= 2;
%             dMin= mean(dBts);
%         end
% 
%         dBts= diag(pdist2(fI2, fNh1));
%         if (mean(dBts) < dMin)
%             face1= 2;
%             face2= 1;
%         end
%     end
%     
%     if (segMainAx == 2)
%         %In this case, by construction of the Bound. Box, we know that the closest faces are 3 and 4 or viceversa
%         fI1= reshape(facets1(3, :, :), [4, 3]);
%         fI2= reshape(facets1(4, :, :), [4, 3]);
% 
%         fNh1= reshape(facets2(3, :, :), [4, 3]);
%         fNh2= reshape(facets2(4, :, :), [4, 3]);
% 
%         dBts= diag(pdist2(fI1, fNh2));
%         if (mean(dBts) < dMin)
%             face1= 3;
%             face2= 4;
%             dMin= mean(dBts);
%         end
% 
%         dBts= diag(pdist2(fI2, fNh1));
%         if (mean(dBts) < dMin)
%             face1= 4;
%             face2= 3;
%         end
%     end
% 
%     if (segMainAx == 3)
%         %In this case, by construction of the Bound. Box, we know that the closest faces are 5 and 6 or viceversa
%         fI1= reshape(facets1(5, :, :), [4, 3]);
%         fI2= reshape(facets1(6, :, :), [4, 3]);
% 
%         fNh1= reshape(facets2(5, :, :), [4, 3]);
%         fNh2= reshape(facets2(6, :, :), [4, 3]);
% 
%         dBts= diag(pdist2(fI1, fNh2));
%         if (mean(dBts) < dMin)
%             face1= 5;
%             face2= 6;
%             dMin= mean(dBts);
%         end
% 
%         dBts= diag(pdist2(fI2, fNh1));
%         if (mean(dBts) < dMin)
%             face1= 6;
%             face2= 5;
%         end
%     end
% end
