%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2016
%Funcion: fGetSegIndices
%Input:
%   allSegments:    structure array containing the segments that were previously computed for all 
%                   the input shapes
%   shapeIx:        index of the Shape for which the graph will be computed
%   shapeName:      (OPTIONAL) name of the Shape for which the graph will be computed
%Output:
%   segIndices:     array (vector) of indices of the Segments in the array
%                   "allSegments" that belong to the same shape, "shapeName"
%
% This function obtains an array containing the indices of the segments 
% (in the array "allSegments"), that were extracted from the same Shape 
% that has the index ShapeIx
%
%--------------------------------------------------------------------------
% function [segIndices] = fGetSegIndicesPerShape(allSegments, shapeIx, shapeName)
function [segIndices] = fGetSegIndicesPerShape(allSegments, shapeIx)
    segIndices= [];
    try        
        for i= 1 : size(allSegments, 2)
            %Find the first segment of the Shape "shapeIx"
%             if (strcmpi(shapeName, allSegments(i).origFileName) == 1)
            if (allSegments(i).ixShape == shapeIx)
                segIndices= [segIndices i];
            end
        end
    catch ME
        segIndices= [];
        errordlg(['Error in fGetSegIndicesPerShape): ' ME.message]);
    end
end