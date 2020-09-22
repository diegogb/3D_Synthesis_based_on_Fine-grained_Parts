%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: May/2016
%Function: fCalcSDFHistoForSegment
%Input:
%   points:     set of 3D points for which the SDF will be computed
%   normals:    array with the normals for the set of 3D points
%   maxSamples: maximum No. of samples that will be used to compute the SDF
%Output:
%   err:        -1 if some error ocurrs; 0, otherwise
%   sdfHisto:   Histogram for the SDF values for the points
% 
% This function estimates the SDF values for a set of points, and obtains a Histogram, 
% capturing the volumetric profile of a segment, based on the method described in Paper (a)
%
% (See also: fGetFarthestPointsFromSeg)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Kaick, Oliver Van, Fish, N., et al. "Shape segmentation by approximate convexity analysis." ACM (TOG) 34.1. (2014)
%   (b) Shapira et al. "Consistent mesh partitioning and skeletonisation using the shape diameter function". (2008)
%(2) This code is based on the implementation developed by NOA FISH for Paper (a). 
%   Code available online at: http://www.cs.tau.ac.il/~noafish/wcseg/
%%%%%
%--------------------------------------------------------------------------
function [err, sdfHisto] = fCalcSDFHistoForSegment(points, normals)
    %Initialize return value with zeros
    err= 0;
    %By default, we use a histogram of size 10
    sdfHisto= zeros(1, 10);
    try
        if (size(points, 1) == size(normals, 1)) && (size(points, 2) == 3) && (size(normals, 2) == 3)
            %We use our own angle value: larger than in Paper (a), but smaller than in Paper (b)
            maxConeAng= (60 / 180) * pi;

            %%%
            %Sample a given number of points from the segment: the number of samples is 
            %computed following Paper (a), inside the function "fGetFarthestPointsFromSeg"
            [err, sampsIx]= fGetFarthestPointsFromSeg(points);
            %%%
            
            if (~err && ~isempty(sampsIx))
                sampPts= points(sampsIx, :);
                
                numSamps = size(sampPts, 1);
                sdfVals = zeros(numSamps, 1);            
                for i = 1 : numSamps
                    dirRays= -normals(sampsIx(i), :);
                    %Calc Rays between the sampled points and all other points of the segment
                    coneRays = bsxfun(@minus, points, sampPts(i, :));
                    %Normalize the rays
                    coneRayNorm = sqrt(sum(coneRays.^2, 2));
                    coneRayNorm = bsxfun(@rdivide, coneRays, coneRayNorm);

                    %Get the angle between the normal of the sampled point and each ray
                    dpRays = sum(bsxfun(@times, coneRayNorm, dirRays), 2);
                    angRays = acos(dpRays);
                    %Find the rays that lie inside the cone
                    raysInCone = find(angRays <= maxConeAng);

                    %If the cone is empty, we try flipping the normals: some of our Input Shapes 
                    %have normals that are in the wrong direction
                    if (isempty(raysInCone))
                        dirRays= -dirRays;
                        dpRays = sum(bsxfun(@times, coneRayNorm, dirRays), 2);
                        angRays = acos(dpRays);

                        raysInCone = find(angRays <= maxConeAng);
                    end

                    if (~isempty(raysInCone))
                        %Get the lengh of the rays that are inside the cone: valid rays
                        valRaysLen= coneRays(raysInCone, :);
                        valRaysLen= sqrt(sum(valRaysLen.^2, 2));
                        %Following Paper (a), weight the length of the valid rays by de cosine
                        %of the angle between the normal of the sample point and the valid rays:
                        %Rays that deviate much will not be considered
                        valRaysLen = valRaysLen .* cos(angRays(raysInCone));
                        [sortedDist, sortedIx]= sort(valRaysLen);

                        chkVsNorms = sum(coneRayNorm(raysInCone,:) .* normals(raysInCone,:), 2);
                        chkVsNorms = chkVsNorms(sortedIx);
                        penalized = find(chkVsNorms <= 0);

                        %Get the median of the valid rays to obtain the SDF                    
                        if (~isempty(penalized))
                            if size(penalized, 1) < size(sortedIx, 1)                            
                                sdfVals(i) = median(sortedDist(setdiff(sortedIx, penalized)));                            
                            end
                        else
                            sdfVals(i) = median(valRaysLen);
                        end
                    end
                end

                %Finally, compute a Histogram of the values of the Sdf
                if (any(sdfVals))
                    sdfHisto= histcounts(sdfVals, 10, 'Normalization', 'probability');            
                end
            end            
        else
            err= -1;
            errordlg('Error in "fCalcSDFHistoForSegment". The size of the required input parameters is not valid');
        end
    catch ME
        err= -1;
        errordlg(['Error computing the SDF for a segment (error in "fCalcSDFHistoForSegment"): ' ME.message]);
    end
end