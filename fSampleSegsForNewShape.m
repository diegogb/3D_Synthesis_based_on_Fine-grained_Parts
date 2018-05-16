%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jul/2017 - XXX/2017
%Funcion: fSampleSegsForNewShape
%Input:
%   segmentsArray:      structure array containing the segments that were previously computed for all 
%                       the input shapes
%   samplingsArray:     structure array containing the samples of the original Shapes
%   segmentsAdj:        structure containing graphs that represent the connectivity 
%                       of the segments for each Shape
%   templateIndices:    indices of the segments of the current Template Shape in the array 
%                       "segmentsArray", (see fGetSegIndicesPerShape)
%   Dist:               Matrix of distances or similarity between each pair of segments, 
%                       (see fGetDistanceBetweenSegs)
%   scoresArray:        array of scores computed for each segment, according to
%                       the distance between adjacent segments on each shape, (see fCalcScoresWithNeighbours)
%   statsPerShape:      array with the mean and std dev. of each shape, 
%                       computed from the values of the "scoresArray". (see fCalcScoresWithNeighbours)
%   minProbThres:       minimum probability threshold for selecting the
%                       replacement of each segment of the Template
%   allowSameShape:     (OPTIONAL) indicates if we sample segments from the
%                       same shape or not(default = true)%   
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%   selSeg:         array with the indices of the segments found as replacements
%                   of the segments of the Template, for synthesizing a new shape
%   selSegScore:    array of scores computed for each replacement Segment
%
% Function that computes the "Energy" formula used to synthesized novel
% shapes, obtaining segments having a high probability of being a good
% replacement for each segment of a Template Shape
%
% (See also fGetDistanceBetweenSegs, fGetSegIndicesPerShape)
%
%%%
% NOTES: 
%   If it is not possible to find any good replacement, the original segment of the template 
%   is kept and returned as the selected segment for synthesizing the new shape
%%%
%--------------------------------------------------------------------------
function [err, selSeg, selSegScore] = fSampleSegsForNewShape(segmentsArray, samplingsArray, segmentsAdj, templateIndices, Dist, ...
    scoresArray, statsPerShape, minProbThres, allowSameShape)
    err= 0;
    selSeg= zeros(size(templateIndices, 2), 1);
    selSegScore= zeros(size(templateIndices, 2), 1);
    try
        if (minProbThres >= 0.1 && minProbThres < 1)
    
            %Compute the first term of the Energy formula: simply, the distance
            %between the segment of the Template and all the other segments
            enerTerm1= zeros(size(templateIndices, 2), size(segmentsArray, 2));
            for i = 1: size(templateIndices, 2)
                enerTerm1(i, :)= Dist(templateIndices(i), :);
            end

            %%%
            %TODO: these thresholds could be parameters of this function, given by the user?
            minSizeThres= 0.6;
            maxSizeThres= 2.0;
            propRatioThres= 1.8;
            
%             allowSameShape= true;
            %Allow or not the replacement of symmetric vs. non-symmetric segments
            allowSymRepNoSym= true;
            %%%
            
            if (~exist('allowSameShape', 'var'))
                allowSameShape= true;
            end
            
            %Get the adjacency graph (neighbours) for the current template
            ixTempl= segmentsArray(templateIndices(1)).ixShape;
            segAdjacency= segmentsAdj(ixTempl).adjacencyMat;
            
            %%%%%
            %We'll use zscores (standarized scores) for making an
            %additional comparison and filtering of the selected segments
            zsTemplate= zscore(scoresArray(templateIndices));
            %%%%%
            
            for i = 1 : size(templateIndices, 2)
                isValidSamp= false;
                
                %Transform the 1st. Term of the Energy into a Probability Density Function
                pT1= exp(-enerTerm1(i, :));
%                 pT1= exp(enerTerm1(i, :));
                pT1= pT1 ./ sum(pT1);
%                 pT1= enerTerm1(i, :) ./ sum(enerTerm1(i, :));
                %Sort the values to keep track of the original indices of the segments
                [~, ixSor]= sort(pT1);
                csPT1= cumsum(pT1);
                %Get all segments with probability higher than the threshold. (If we have a segment with probability= 1, 
                %that must be the same current segment from the template, so we exclude that one)
                allPosSamps= find(csPT1 > minProbThres & csPT1 < 1.0);
                
                %%%%%%%%%                
                %Do not uncomment these lines!!: they're only for
                %experiments to show what is the energy value for a "bad case"
%                 minProbThres= 0.05;
%                 allPosSamps= find(csPT1 > minProbThres & csPT1 < 0.5);
                %%%%%%%%%
                
                if length(allPosSamps) > 1                
                    %Get a random permutation of length equal to the number
                    %of segments with probability higher than the threshold
                    ixrs= randperm(length(allPosSamps));
                    %Now, start with the first value of the permutation and check if its a valid segment, 
                    %according to our additional criteria below
                    sampIx= allPosSamps(ixrs(1));
                else
                    %In case we don't have any segment with probability
                    %higher than the threshold: this should not happen
                    sampIx= ixSor(templateIndices(i));
                    isValidSamp= true;
                end

                ptsRef= samplingsArray(segmentsArray(templateIndices(i)).ixShape).samples(segmentsArray(templateIndices(i)).ixSamples, :);
                szRef= max(ptsRef) - min(ptsRef);
                
                %Check if sampled segment is valid; otherwise, keep looking for a possible 
                %replacement until a given number of iterations have been executed
                numT= 1;
                while (~err && ~isValidSamp && numT < length(allPosSamps))                    
                    %Verify if selected segment is not self:
                    if (ixSor(sampIx) ~= templateIndices(i))
                        ptsNew= samplingsArray(segmentsArray(ixSor(sampIx)).ixShape).samples(segmentsArray(ixSor(sampIx)).ixSamples, :);
                        %Check the size threshold
                        szNew= max(ptsNew) - min(ptsNew);

                        maxRatP1= max(max(abs((szRef(1) / szRef(2))-(szNew(1) / szNew(2))), abs((szRef(1) / szRef(3))-(szNew(1) / szNew(3)))), ...
                            abs((szRef(2) / szRef(3))-(szNew(2) / szNew(3))));
                        maxSCFact= szRef(segmentsArray(templateIndices(i)).segmentMainAx) / szNew(segmentsArray(ixSor(sampIx)).segmentMainAx);
    
                        if (maxRatP1 <= propRatioThres) && (maxSCFact >= minSizeThres && maxSCFact <= maxSizeThres)
                            %Check if selected segment is not from same symmetry group:
                            [err, ~, symGroupIndices] = fGetSymmetryChain(templateIndices(i), segmentsArray, templateIndices);
                            if (~err && ~isempty(symGroupIndices))
                                if (~err && ~ismember(ixSor(sampIx), symGroupIndices))
                                    %Check if reference segment vs. new segment is symmetric                
                                    if (allowSymRepNoSym || ~isempty(segmentsArray(ixSor(sampIx)).symLabels))
                                        %Check if segments is from the same shape
                                        if (allowSameShape || segmentsArray(ixSor(sampIx)).ixShape ~= segmentsArray(templateIndices(i)).ixShape)
%                                             isValidSamp= true;
                                            
                                            newShp= segmentsArray(ixSor(sampIx)).ixShape;
                                            zscoreOfSamp= (scoresArray(ixSor(sampIx)) - statsPerShape(newShp).mean) / statsPerShape(newShp).stdd;
                                          
                                            if (abs(zscoreOfSamp - zsTemplate(i)) <= 2)
                                                %We filter out the segment if it's more than 2 standard deviations apart
                                                %in comparison to the segment of the template
                                                isValidSamp= true;
                                            end
                                            %%%
                                        end
                                    end
                                end
                            elseif (~err)
                                %Check if reference segment vs. new segment is symmetric                
                                if (allowSymRepNoSym || isempty(segmentsArray(ixSor(sampIx)).symLabels))
                                    %Check if segments is from the same shape
                                    if (allowSameShape || segmentsArray(ixSor(sampIx)).ixShape ~= segmentsArray(templateIndices(i)).ixShape)
%                                         isValidSamp= true;
                                        
                                        newShp= segmentsArray(ixSor(sampIx)).ixShape;
                                        zscoreOfSamp= (scoresArray(ixSor(sampIx)) - statsPerShape(newShp).mean) / statsPerShape(newShp).stdd;

                                        if (abs(zscoreOfSamp - zsTemplate(i)) <= 2)
                                            %We filter out the segment if it's more than 2 standard deviations apart
                                            %in comparison to the segment of the template
                                            isValidSamp= true;
                                        end
                                    end
                                end                
                            else
                                break;
                            end
                        end
                    end

                    if (~isValidSamp)
                        %Get the next candidate sample from the permutation
                        numT= numT + 1;
                        if (numT <= length(allPosSamps))
                            sampIx= allPosSamps(ixrs(numT));
                        else
                            %%%
                            %If we don't find any segment with enough similarity, we keep the same segment of the Template
                            % THIS SHOULD NOT HAPPEN
                            sampIx= ixSor(templateIndices(i));
                            isValidSamp= true;
                            %%%
                        end
                    end
                end
                %%%

                %Get the index of the segment corresponding to the sampled value: this is the segment 
                %that will be used as replacement of the segment of the Template
                selSeg(i)= ixSor(sampIx);
                %Save the probability associated to the replacement, as initial
                %score for the new segment
                selSegScore(i)= csPT1(sampIx);            
            end

            %Compute the 2nd. Term of the Energy formula...
            enerTerm2= zeros(size(segmentsArray, 2), size(segmentsArray, 2));
            %For each segment "i" of the template shape...
            for i = 1 : size(templateIndices, 2)
                ixNeighs= find(segAdjacency(i, :) == 1);

                for nh= 1 : length(ixNeighs)
                    %We get the distance between the segment "i" of the template and each one of its neighbors
                    enerBetweenNeighs= Dist(templateIndices(i), templateIndices(ixNeighs(nh)));

                    %We get now the distance between the segment that will replace the segment "i" 
                    %and the segment that will replace its neighbor
                    distNewNeighs= Dist(selSeg(i), selSeg(ixNeighs(nh)));
                    
                    %Finally, we get the difference between the two distances obtained above
                    enerTerm2(selSeg(i), selSeg(ixNeighs(nh)))= abs(enerBetweenNeighs - distNewNeighs);
                end
            end

            %Compute the PDF for the 2nd. energy term...
            pT2= exp(-enerTerm2);
            pT2= pT2 / size(segmentsArray, 2);            
            csPT2= cumsum(pT2);

            %Sum up the probability of the adjacent segments
            for i = 1 : size(templateIndices, 2)    
                ixNewNeighs= selSeg(segAdjacency(i, :) == 1);
                if size(ixNewNeighs, 1)> 0
                    %TODO: check if this is correct!!!:
                    %We divide by the number of neighbors to get a value normalized in the range [0, 1]
                    selSegScore(i)= selSegScore(i) + (sum(csPT2(selSeg(i), ixNewNeighs)) / size(ixNewNeighs, 1));
                end
            end

            %Finally, make sure that we preserve the symmetries according to
            %the Template Shape
            for i = 1 : size(templateIndices, 2)
                [err, symmLbl] = fGetSymmetryChain(templateIndices(i), segmentsArray, templateIndices);
                if (~isempty(symmLbl))
                    %Sort the symmetric replacement segments by score (probability), and select as
                    %replacement for all the symmetric segments, the one having the maximum score
                    [maxSC, maxIx]= max(selSegScore(symmLbl));
                    selSeg(symmLbl)= repmat(selSeg(symmLbl(maxIx)), length(maxIx), 1);
                    selSegScore(symmLbl)= repmat(maxSC, length(maxIx), 1);
                end
            end 
        else
            sErr= 'Error sampling segments to synthize a novel shape: the probability threshold must be a valid value in the range [0.5, 1)';
            errordlg(sErr);
        end
    catch ME
        err= -1;
        sErr= ['Error sampling segments to synthize a novel shape (error in fSampleSegsForNewShape): ' ME.message];
        errordlg(sErr);
    end
end