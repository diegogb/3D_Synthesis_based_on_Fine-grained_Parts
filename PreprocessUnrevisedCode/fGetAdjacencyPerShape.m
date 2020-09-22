%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jul/2017 - XXX/2017
%Function: fGetAdjacencyPerShape
%Input:
%   segmentsArray:  structure array containing the segments that were previously computed for all 
%                   the input shapes
%   samplingsArray: array of structures containing the sampled points from the Meshes
%Output:
%   segmentsAdj:    structure array with one element per shape, containing a graph containing the
%                   connectivity between segments, and information about
%                   the main axis of orientation of each segment
%
% This function computes graphs representing the connectivity between
% segments for each Shape contained in the array "samplingsArray"
%
% (See also fGetGraphForSegments)
%
%--------------------------------------------------------------------------
function [err, segmentsAdj] = fGetAdjacencyPerShape(segmentsArray, samplingsArray)
    %Initialize output
    segmentsAdj= struct([]);
    try
        numShapes= size(samplingsArray, 2);
        for i= 1 : numShapes
%             [err, adjMat, adjAddInfo] = fGetGraphForSegments(segmentsArray, samplingsArray, i, segSymArray);
            [err, adjMat, adjAddInfo] = fGetGraphForSegments(segmentsArray, samplingsArray, i);
%             [err, adjMat, adjAddInfo, segMainAx] = fGetGraphForSegments(segmentsArray, samplingsArray, i);
            if (~err)                
                adj.adjacencyMat= adjMat;
                adj.adjacencyInfo= adjAddInfo;
%                 adj.segmentsMainAx= segMainAx;
%                 adj.closestPoints= closest;
                segmentsAdj= [segmentsAdj adj];
            else
                segmentsAdj= struct([]);
                break;
            end            
        end
    catch ME
        err= -1;
        segmentsAdj= struct([]);
        errordlg(['Error computing the graph for the segments of the input shapes (error in "fGetAdjacencyPerShape"): ' ME.message]);
    end    
end