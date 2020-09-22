%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2017 - XXX/2017
%Function: fGetSubDividedMeshes
%Input:
%   meshesOrig:     array of structures with the original meshes
%   maxEdgePercent: (OPTIONAL) value between 0 and 1 representing the percentage of
%                   the maximum edge length that is going to be used as
%                   threshold for subdividing the meshes. Default= 0.05
%   meshNorms:      (OPTIONAL) normals of the original meshes, which can be
%                   transfered (copied) to the faces of the subdivided mesh
%Output:
%   meshesArraySD:  array of structures with the subdivided meshes
%
% This function returns an array of structures containing one subdivided mesh for each mesh in the array meshesOrig. 
% The subidivision is performed until no face has any edge longer than the threshold "maxEdgePercent"
%%%%%
% NOTE: this function is basically a wrapper for the function
% "mesh_subdivide", developed by Professor Oliver van Kaick
%%%%%
%
% (See also mesh_subdivide)
%--------------------------------------------------------------------------
function [err, meshesArraySD] = fGetSubDividedMeshes(meshesOrig, maxEdgePercen, meshNorms)
    %Initialize output
    meshesArraySD= struct([]);
    try
        if (~exist('maxEdgePercen', 'var'))
            maxEdgePercen = 0.05;
        else
            if (maxEdgePercen < 0) || (maxEdgePercen > 1)
                maxEdgePercen= 0.05;
            end
        end
        
        numShapes= size(meshesOrig, 2);
        for i= 1 : numShapes
            M= meshesOrig(i);
            allLens= zeros(size(M.faces, 1), 3);
            
            %Get the length of the edges of the mesh
            e1= M.vertices(M.faces(:, 1), :) - M.vertices(M.faces(:, 2), :);
            lenE1= sqrt(sum(e1 .^ 2, 2));
            e2= M.vertices(M.faces(:, 2), :) - M.vertices(M.faces(:, 3), :);
            lenE2= sqrt(sum(e2 .^ 2, 2));
            e3= M.vertices(M.faces(:, 3), :) - M.vertices(M.faces(:, 1), :);
            lenE3= sqrt(sum(e3 .^ 2, 2));
            
            allLens(:, 1)= lenE1;
            allLens(:, 2)= lenE2;
            allLens(:, 3)= lenE3;
            longestEdge= max(max(allLens));
            shortestEdge= min(min(allLens));
            
            minLen= (longestEdge - shortestEdge) * maxEdgePercen;
            
            if (exist('meshNorms', 'var'))
                mNormals= meshNorms(i).normals;
                Msub= mesh_subdivide(M, minLen, mNormals);
            else
                Msub= mesh_subdivide(M, minLen);
            end
            
            meshesArraySD = [meshesArraySD Msub];
        end
    catch ME
        err= -1;
        meshesArraySD= struct([]);
        errordlg(['Error computing subdivided meshes (error in "fGetSubDividedMeshes"): ' ME.message]);
    end    
end