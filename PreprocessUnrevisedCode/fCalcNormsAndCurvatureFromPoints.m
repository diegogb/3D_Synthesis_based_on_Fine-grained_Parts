%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jun/2017
%Function: fCalcNormsAndCurvatureFromPoints
%Input:
%   points: set of 3D points for which the normals will be obtained
%   k:      (OPTIONAL) number of neighbors that will be used in the
%           calculation (Default = 6)
%Output:
%   normals:            Array with the computed normals for the set of 3D points
%   surfaceVariation:   returns the estimated surface variation at each point of the input set of points
%
% This function approximates the normals of a set of 3D points, by using a set of k neighbors for each point 
% and performing PCA over the covariance matrix computed for the points of this neighborhood; see References (a), (b).
% Also, the eigenvalues obtained by the PCA are used to compute the surface variation, following Reference (b)
%
% (See also fGetDistanceBetweenSegs)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Hoppe, H., et al. "Surface reconstruction from unorganized points". Vol. 26. No. 2. ACM, (1992)
%   (b) Pauly, M., et al. "Efficient simplification of point-sampled surfaces." Proc. of the Conf. on Visualization'02. IEEE Comp. Soc, (2002).
%   (c) http://pointclouds.org/documentation/tutorials/normal_estimation.php
%%%%%
%--------------------------------------------------------------------------
function [normals, surfaceVariation] = fCalcNormsAndCurvatureFromPoints(points, neighSize)
    %Initialize return value with zeros
    normals= zeros(size(points, 1), 3);
%     curvature= zeros(size(points, 1), 1);
    surfaceVariation= zeros(size(points, 1), 1);
    try
        if (size(points, 1) > 3) && (size(points, 2) == 3)
            numPoints= size(points, 1);
            if (~exist('neighSize', 'var'))
                neighSize= 6;
            else
                if (neighSize < 6)
                    neighSize= 6;
                elseif (neighSize > 40) || (neighSize > size(points, 1))
                    neighSize= min(40, size(points, 1));
                end
            end
            
            %Define a kdtree to search the "k" nearest neighbours for each point
            kdt= KDTreeSearcher(points);            
            neighIdx= knnsearch(kdt, points, 'k', neighSize + 1);
            %Function knnsearch returns the self point as first closest point, so remove it
            neighIdx(:, 1)= [];
            
            for i= 1 : numPoints
                %For each point get the covariance matrix for the nearest neighbours
                matCCov = cov(points(neighIdx(i, :), :));
                
                %Get eigenvalues and eigenvectors for the covariance matrix
                [eVecs, eVals]= eigs(matCCov);
                eVals= diag(eVals);
                
                %According to References (a), (b), (c), the eigenvector associated 
                %to the min. eigenvalue is the approximation to the normal
                %for the patch given by the nearest neighbours
                [minEv, minIx] = min(eVals);
                normals(i, :)= eVecs(:, minIx);
                %Following Reference (b), (c):
                surfaceVariation(i)= minEv / sum(eVals);
            end
            normals = round(normals, 4);
            
%             %Instead of the procedure presented in Reference (a), we try to find normals 
%             %consistenly oriented, just by using the "most common" normal of the neighbors
%             for i= 1 : numPoints
%                 [uniN, ~, ixUN] = unique(normals(neighIdx(i, :), :), 'rows');
%                 
%                 [~, ixMaxN] = max(accumarray(ixUN, 1));
%                 normals(i, :)= uniN(ixMaxN, :);
%             end
        else
            errordlg('Error in "fCalcNormsAndCurvatureFromPoints". The size of the required input parameter is not valid');
        end        
    catch ME
        errordlg(['Error computing the overall geometry for a set of points (error in "fCalcNormsAndCurvatureFromPoints"): ' ME.message]);
    end
end