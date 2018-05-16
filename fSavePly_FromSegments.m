%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Feb/2016-Abr/2016
%Function: fSavePly_FromSegments
%Input:
%   patchesArray:   Array containing the points that will be saved in a single PLY file   
%   pathToSave:     name (including path), of the file where the point sampling will be saved
%   colorsPerPatch: (optional) 3D vertex that represents a color that will be associated to all the points of a given segment
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%   
% This function goes trough all the elements of an Array that contains
% several structures with "point clouds" that are Segments taken from a
% larger "entire" shape (or from a set of Shapes)
%--------------------------------------------------------------------------
function err= fSavePly_FromSegments(patchesArray, pathToSave, colorsPerPatch)
    allSegments = [];
    allNorms= [];
    allColors= [];
    try
        %For each segment...
        for i= 1 : length(patchesArray)
%             allSegments = cat(1, allSegments, patchesArray(i).seg);
            allSegments = cat(1, allSegments, patchesArray(i).points);
            %cat(3, IMAGES, Img);
            allNorms= cat(1, allNorms, patchesArray(i).normals);
            
            if (nargin == 3)
                %If there are colors defined per each patch (segment),
                %"repeat" the color for each point in the segment                
                if (exist('colorsPerPatch', 'var'))
%                     allColors= cat(1, allColors, repmat(colorsPerPatch(i, :), size(patchesArray(i).seg, 1), 1));
                    allColors= cat(1, allColors, repmat(colorsPerPatch(i, :), size(patchesArray(i).points, 1), 1));
                end
            end
        end
        
        if (~isempty(allColors))
            err = points_write_ply(pathToSave, allSegments, allNorms, allColors);
        else
            err = points_write_ply(pathToSave, allSegments, allNorms);
            %err = points_write_ply(pathToSave, allSegments);
        end
        
    catch ME
        err= -1;
        sErr= ['Error writing the "PLY" file using the selected fine-grained segments. (error in fSavePly_FromSegments)' ME.message];
        errordlg(sErr);
    end
end