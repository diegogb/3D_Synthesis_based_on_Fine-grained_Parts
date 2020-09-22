%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: May/2016
%Function: fCalcSDF
%Input:
%   theMesh:        Structure with the faces and vertices of the Mesh
%   bariCent:       Matrix with the baricentric Coords of each Face of the Mesh.
%                   This points will be used as origins of the rays that will be used to compute the volume
%   normals:        Matrix with the normals of each Face of the Mesh. The rays used to compute the volume,
%                   will point towards the opposite side of this vector
%   doNormalize:    (OPTIONAL) Boolean value indicating if the values obtained should be normalized, 
%                   where the normalization is performed as indicated in Reference (a)
%Output:
%   sdfVals:    Array with the computed SDF value of each Face of the Mesh
% 
% This function computes the SDF values of a Mesh, i.e: the approx. volume at each face of a Mesh
%(See also fCalcIntersectRayMesh, fCalcRandPointsOnCone)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Shapira et al. "Consistent mesh partitioning and skeletonisation using the shape diameter function". (2008)
%(2) This function computes the SDF per face, and does not include a smoothing as the original method from (a)
%%%%%
%--------------------------------------------------------------------------
function [sdfVals] = fCalcSDF(theMesh, baris, normals, areas, doNormalize)
    %Initialize return value with zeros
    sdfVals = zeros(size(theMesh.faces, 1), 1);
    
    zeroEpsil= 1e-4;
    %Number of rays per face, (rays contained inside the cone): value as suggested in paper (a)
    numRaysPerFace= 30;    
    try
        %Compute vectors representing two edges of each triangle in the Mesh. 
        %This is used later for computing the Ray-Mesh intersections
        e1= theMesh.vertices(theMesh.faces(:, 2), :) - theMesh.vertices(theMesh.faces(:, 1), :);
        e2= theMesh.vertices(theMesh.faces(:, 3), :) - theMesh.vertices(theMesh.faces(:, 1), :);
        
        for i= 1 : size(theMesh.faces, 1)
%         parfor i= 1 : size(theMesh.faces, 1)
            tempVals= zeros(numRaysPerFace + 1, 1);
            if areas(i) > zeroEpsil                
                %For each face, start with a ray exactly towards the opposite direction of the normal
                [errI, numHitsI, ixHit, isecLen]= fCalcIntersectRayMesh(theMesh, e1, e2, baris(i, :), -normals(i, :), i, true);
                if (~errI)
                    if numHitsI > 0
                        tempVals(1)= isecLen(ixHit);
                    end
                end

                %For the cone of rays: get the rays towards the opposite direction of the normal
    %             coneRayDirs= fCalcRandPointsOnCone(numRaysPerFace, maxConeAng, baris(i, :), -normals(i, :));
%                 coneRayDirs= fCalcRandPointsOnCone(numRaysPerFace, -normals(i, :));
                coneRayDirs= fCalcRandPointsOnCone(-normals(i, :));

                %Cast rays from the baricenter of the face, towards directions inside the cone
                %and calc the first intersection of each ray with the Mesh
                for j= 1 : numRaysPerFace
                    [errI, numHitsI, ixHit, isecLen]= fCalcIntersectRayMesh(theMesh, e1, e2, baris(i, :), coneRayDirs(j, :), i, true);
                    if (~errI)
                        if numHitsI > 0
                            tempVals(j + 1)= isecLen(ixHit);
                        end
                    end
                end

%                 %Compute average of interesection lengths
%                 meanLen= mean(tempVals);
                
                %Compute MEDIAN of interesection lengths
                medianLen= median(tempVals);
                
                %Look for rays for which the intersection distance is too large in comparison with the mean length, and discard them
                stDev= std(tempVals);
%                 validVals= tempVals(abs(tempVals - meanLen) <= 0.5 * stDev);
                validVals= tempVals(abs(tempVals - medianLen) <= 0.5 * stDev);
                if (size(validVals, 1)) > 0 
                    sdfVals(i) = sum(validVals) / size(validVals, 1);
                end
                
                %TODO: Include "outliers" removal (??)
            end
        end
        %For values too close to zero, assign zero
        sdfVals(sdfVals >= -zeroEpsil & sdfVals <= zeroEpsil)= 0;
        
        if exist('doNormalize', 'var')
            doNormalizeVals= doNormalize;
        else
            doNormalizeVals= false;
        end
        
        if doNormalizeVals
            %Normalization factor
            normAlpha= 4;
            sdfVals= log(((sdfVals - min(sdfVals)) / (max(sdfVals) - min(sdfVals))) * normAlpha + 1) / log(normAlpha + 1);
        end
    catch ME
        errordlg(['Error computing the SDF values of the Mesh (error in fCalcSDF): ' ME.message]);
    end
end
