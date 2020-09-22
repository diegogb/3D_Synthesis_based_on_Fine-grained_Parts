%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Aug/2016 - XXX/2017
%Function: fCalcOverallGeomFromEigen
%Input:
%   points:	set of 3D points for which a vector representing the overall geometry 
%           of the points will be computed, using PCA (eigenvalues): see Reference (a)
%Output:
%   eigenGeomVector:    Array with the computed vector representing the
%                       overall geometry of the segment, according to Paper (a)
%   eVecs:              matrix containing in its columns the 3 principal eigenvectors for the segment
%   eVals:              3-component vector (vector-column) containing the first 3 eigenvalues for the segment
%
% This function obtains a 3D vector that represents the overall geometry of a set of points, 
% by means of computing the three principal eigenvalues for the matrix of points
%
% (See also fGetDistanceBetweenSegs)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Sidi, O., van Kaick, O., et al. "Unsupervised co-segmentation of a set of shapes via descriptor-space spectral clustering". (2011)
%%%%%
%--------------------------------------------------------------------------
function [eigenGeomVector, eVecs, eVals] = fCalcOverallGeomFromEigen(points)
    %Initialize return value with zeros
    eigenGeomVector = zeros(1, 3);
    try
        if (size(points, 1) > 3) && (size(points, 2) == 3)
            matCCov = cov(points);
            [eVecs, eVals]= eigs(matCCov);
            eVals= diag(eVals);
            if size(eVals, 1) == 3
                %Sort eigenValues and eigenVectors
                [eVals, ixev] = sort(eVals, 'descend');
                eVecs= eVecs(:, ixev);
                
                %Using the eigenvalues, compute each of the components of the
                %"overall geometry" vector, according to Paper (a)
                mu_l= (eVals(1) - eVals(2)) / sum(eVals);
                mu_p= (2 * (eVals(2) - eVals(3))) / sum(eVals);
                mu_s= (3 * eVals(3)) / sum(eVals);

                eigenGeomVector= [mu_l mu_p mu_s];
            else
                errordlg('Error in "fCalcOverallGeomFromEigen": there was an error performing the eigendecomposition');
            end
        else
            errordlg('Error in "fCalcOverallGeomFromEigen". The size of the required input parameter is incorrect: at least three 3D points are required');
        end        
    catch ME
        eigenGeomVector = zeros(1, 3);
        errordlg(['Error computing the overall geometry for a set of points (error in "fCalcOverallGeomFromEigen"): ' ME.message]);
    end    
end