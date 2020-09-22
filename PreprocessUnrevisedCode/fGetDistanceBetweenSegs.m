%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jun/2016-XXX/2017
%Funcion: fGetDistanceBetweenSegs
%Input:
%   allSegments:    array containing the segments to compare, i.e. the segments 
%                   for which the distance or similarity will be computed
%   samplingsArray: array containing the points that were sampled from the input 3D Shapes, 
%                   (i.e., the Point Clouds from where the segments were obtained)
%   doNormalize:    normalize the final distance values in the range 0-1
%Output:
%   err:        -1 if some error ocurrs; 0, otherwise
%   Dist:       Pairwise Distance matrix between each pair of segments
%
% This function compares the similarity or distance between each pair of
% parts or segments given by the array "allSegments". 
%
% (See also fCalcAnglesForPFH, fCalcOverallGeomFromEigen, fCalcSDFHistoForSegment, fCalcNormsAndCurvatureFromPoints, fCalcPointDistriHistoForSegment)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Sidi, O., van Kaick, O., et al. "Unsupervised co-segmentation of a set of shapes via descriptor-space spectral clustering". (2011)
%   (b) https://en.wikipedia.org/wiki/Earth_mover%27s_distance
%   (c) Lun, Z., et al. "Functionality preserving shape style transfer." ACM (TOG) 35.6 (2016).
%%%%%
%--------------------------------------------------------------------------
function [err, Dist] = fGetDistanceBetweenSegs(allSegments, samplingsArray, doNormalize)
    try
        numSegs= size(allSegments, 2);
        numShapes= size(samplingsArray, 2);
        
        heightsPerShape= zeros(size(samplingsArray(1).samples, 1), numShapes);
        %%%
        centsPerShape= zeros(numShapes, 3);
        %%%
        
        hCurv = zeros(numSegs, 10);
        hSDF = zeros(numSegs, 10);
        hNormals= zeros(numSegs, 125);
        eigGV = zeros(numSegs, 3);
        heightMS = zeros(numSegs, 1);
        minP = zeros(numSegs, 3);
        maxP = zeros(numSegs, 3);
        centMS = zeros(numSegs, 1);
        hPDist= zeros(numSegs, 64);
        %%%
        
        %%%        
%         valNeighConf= zeros(numSegs, 1);
        %%%
        
        Dist = zeros(numSegs);
        
        %For each point, get the height from the base of the Shape
        for i = 1 : numShapes
            %Height from the base per shape:
            heightsPerShape(:, i)= fCalcHeightfromBase(samplingsArray(i).samples);
            %Centroid of each shape:
            centsPerShape(i, :)= mean(samplingsArray(i).samples);
        end
 
%         tic;
        for i = 1 : numSegs
            segPoints = samplingsArray(allSegments(i).ixShape).samples(allSegments(i).ixSamples, :);
            segNormals = samplingsArray(allSegments(i).ixShape).normals(allSegments(i).ixSamples, :);

            %Get the average Height from the Base for each segment
            heightMS(i) = mean(heightsPerShape(allSegments(i).ixSamples, allSegments(i).ixShape));

            %New descriptors, based on Paper (c): highest and lowest point per segment
%             minH(i)= min(heightsPerShape(allSegments(i).ixSamples, allSegments(i).ixShape));
%             maxH(i)= min(heightsPerShape(allSegments(i).ixSamples, allSegments(i).ixShape));
            minP(i, :)= min(segPoints);
            maxP(i, :)= max(segPoints);
            
            %Get the centroid of each segment
            meanSeg= mean(segPoints);
            %Our version of a descriptor from Paper (c): distance between projection of centroid of each shape to the ground plane
            %and projection of centroid of each segment to the ground plane: assuming the ground plane is always on Y= -1
            centMS(i)= sum(([meanSeg(1) -1 meanSeg(3)] - [centsPerShape(allSegments(i).ixShape, 1) -1 centsPerShape(allSegments(i).ixShape, 3)]) .^ 2, 2);

            %Compute a histogram of 10 bins for the estimated Gaussian Curvature
            [~, curvature] = fCalcNormsAndCurvatureFromPoints(segPoints, 10);
            hCurv(i, :)= histcounts(curvature, 10, 'Normalization', 'probability');
            
            %%%
            [err, hSDF(i, :)]= fCalcSDFHistoForSegment(segPoints, segNormals);
            if (~err)
                [err, hNormals(i, :)] = fCalcAnglesForPFH(segPoints, segNormals);
                if (~err)
                    eigGV(i, :)= fCalcOverallGeomFromEigen(segPoints);
                    if (any(eigGV(i, :)))
                        %Based on Paper (c): compute a histogram for the estimated distribution 
                        %of points (mass distribution) of the segment
                        [err, hPDist(i, :)]= fCalcPointDistriHistoForSegment(segPoints);
                        if (err)
                            break;
                        end
                    else
                        err= -1;
                        break;
                    end
                else
                    break;
                end
            else
                break;
            end
%             end
        end
%         toc;
        
        if (~err)
%             wGC= 0.2;
%             wCurv= 0.125;
%             wSDF= 0.125;
%             wNorm= 0.125;
%             wEigs= 0.125;
%             wHeight= 0.125;
%             %%%
%             wAngPCA = 0.125;
%             wCent = 0.125;
%             %%%
%             wSize= 0.125;
%             wCurv= 1;
%             wSDF= 1;
%             wNorm= 1;
%             wEigs= 1;
%             wHeight= 1;
%             %%%
%             %For new descriptors based on Paper (c)
%             wAngPCA = 1;
% %             wCent = 1;
%             wPD= 1;
%             wCentRelS= 1;
            
%             tic;
            for i = 1 : numSegs
                for j= i + 1 : numSegs
                    dTermCurv = earthMoverDist(hCurv(i, :), hCurv(j, :)) ^ 2;
                    dTermSDF = earthMoverDist(hSDF(i, :), hSDF(j, :)) ^ 2;
                    dTermNormals = earthMoverDist(hNormals(i, :), hNormals(j, :)) ^ 2;

                    dSegEigs= sum((eigGV(i, :) - eigGV(j, :)) .^ 2, 2);
                    dSegHeight= abs(heightMS(i) - heightMS(j)) ^ 2;
                    
                    %%%
                    dSegMin= sum((minP(i, :) - minP(j, :)) .^ 2, 2);
                    dSegMax= sum((maxP(i, :) - maxP(j, :)) .^ 2, 2);
                    %%%
                    
                    %%%
                    dSegCentRS= abs(centMS(i) - centMS(j)) ^ 2;
                    dTermPD = earthMoverDist(hPDist(i, :), hPDist(j, :)) ^ 2;
                    %%%

                    %%%%%
                    %We do not consider the "Angle between the principal axis, since it seems to be inconsistent"
%                     Dist(i, j) = sqrt(dTermCurv + dTermSDF + dTermNormals + dSegEigs + dSegHeight + dSegMinH + dSegMaxH + dSegCentRS + dSegAng + dTermPD);

                    Dist(i, j) = sqrt(dTermCurv + dTermSDF + dTermNormals + dTermPD + dSegEigs + dSegCentRS + dSegHeight + dSegMin + dSegMax);
                    
%                     Dist(i, j) = sqrt(dTermCurv + dTermSDF + dTermNormals + dTermPD + dSegEigs + dSegCentRS + dSegHeight + dSegMin + dSegMax + dConfNeigh);
                    %%%%%
                                        
                    Dist(j, i)= Dist(i, j);
                end
            end
%             toc;
            
            %Normalize the distance between 0 and 1
            if (exist('doNormalize', 'var'))
                if (doNormalize)
                    Dist= Dist / max(max(Dist));
                end
            end
        end
        
    catch ME
        err= -1;
        sErr= ['Error computing the similarity or distance between the segments (error in fGetDistanceBetweenSegs): ' ME.message];
        errordlg(sErr);
    end
end