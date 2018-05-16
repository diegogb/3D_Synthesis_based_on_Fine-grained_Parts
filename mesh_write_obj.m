%
% This function writes a mesh structure to an obj file
%
% function mesh_write_obj(M, filename);
%
% Input -
%   - M: triangle mesh: M.vertices(i, :) represents the 3D coordinates
%   of vertex 'i', while M.faces(i, :) contains the indices of the three
%   vertices that compose face 'i'. If the mesh also has color
%   attributes defined in M.FaceVertexCData, they are written to the
%   file as entries of the form 'vc <r> <g> <b>' or 'fc <r> <g> <b>'.
%   - filename: name of obj file to create
%
% Output -
%   - status: this variable is 0 if the file was succesfuly written, or
%   1 otherwise
%
%   The output mesh is similar to that written by the mesh_write_smf
%   function, but the way of specifying colors changes.
%
% See also mesh_write, mesh_write_smf
%
function status = mesh_write_obj(M, filename)
%
% Copyright (c) 2008, 2009, 2010 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Open file
fid = fopen(filename, 'w');
status = 0;
if fid == -1
    disp(['ERROR: could not open file "' filename '"']);
    status = 1;
    return;
end

% Write vertices
for i = 1:size(M.vertices, 1)
    fprintf(fid, 'v %f %f %f\n', ...
                M.vertices(i, 1), M.vertices(i, 2), M.vertices(i, 3));
end

% Write faces
for i = 1:size(M.faces, 1)
    fprintf(fid, 'f %d %d %d\n', ...
                M.faces(i, 1), M.faces(i, 2), M.faces(i, 3));
end

% Write colors, if defined
if (isfield(M, 'FaceVertexCData'))
    % Check if colors are defined as RGB or a function
    if size(M.FaceVertexCData, 2) == 1
        % Map function values to RGB data
        % Map function values from [min, max] to [0, 2/3], which is red
        % to blue hue in an HSV color map
        func = mesh_map_val(M.FaceVertexCData, 0, 2/3);
        % Use max(func)-func as the hue value and get RGB colors
        % Saturation and brightness are set to 0.8
        clr = hsv2rgb([max(func)-func ones(length(func), 1)*0.8 ...
                            ones(length(func),1)*0.8]);
    else
        % If we have RGB data, just use it
        clr = M.FaceVertexCData;
    end
    % Convert colors to uchar
    for i = 1:size(clr, 1)
        for j = 1:3
            clr(i, j) = floor(clr(i,j)*255);
        end
    end
    % Print colors
    if size(clr, 1) == size(M.vertices, 1)
        for i = 1:size(M.FaceVertexCData, 1)
            fprintf(fid, 'vc %d %d %d\n', ...
                    clr(i, 1), clr(i, 2), clr(i, 3));
        end
    else
        for i = 1:size(M.FaceVertexCData, 1)
            fprintf(fid, 'fc %d %d %d\n', ...
                    clr(i, 1), clr(i, 2), clr(i, 3));
        end
    end
end

% Close file
fclose(fid);
