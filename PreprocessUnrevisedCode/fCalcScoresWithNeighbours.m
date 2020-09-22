%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jul/2017 - XXX/2018
%Funcion: fCalcScoresWithNeighbours
%Input:
%   segmentsArray:      structure array containing the segments that were previously computed for all 
%                       the input shapes
%   segmentsAdj:        structure containing graphs that represent the connectivity 
%                       of the segments for each Shape
%   Dist:               Matrix of distances or similarity between each pair of segments, 
%                       (see fGetDistanceBetweenSegs)
%   
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%   statsPerShape:  array with the mean and std dev. of each shape, computed from the values of the "scoresArray"
%   scoresArray:    array of scores computed for each segment
%
% Function that computes a total score as the distance or similarity value between all the neighbours of each segment,
% as well as the mean and standard deviation of these scores for each shape
%
% (See also fGetDistanceBetweenSegs, fGetSegIndicesPerShape, fSampleSegsForNewShape)
%
%--------------------------------------------------------------------------
function [err, statsPerShape, scoresArray] = fCalcScoresWithNeighbours(segmentsArray, segmentsAdj, Dist)
    err= 0;    
    try
        numShapes= size(segmentsAdj, 2);
        statsPerShape= struct([]);
        scoresArray= zeros(size(segmentsArray, 2), 1);
        for ns= 1 : numShapes
            segsOfShape= fGetSegIndicesPerShape(segmentsArray, ns);
            segAdjacency= segmentsAdj(ns).adjacencyMat;
            
            for i = 1 : size(segsOfShape, 2)
                ixNeighs= find(segAdjacency(i, :) == 1);
                scoresArray(segsOfShape(i))= sum(Dist(segsOfShape(i), segsOfShape(ixNeighs)));
            end
            
            statsPerShape(ns).mean= mean(scoresArray(segsOfShape));
            statsPerShape(ns).stdd= std(scoresArray(segsOfShape));
        end
    catch ME
        err= -1;
        sErr= ['Error in fCalcDistValPerShape: ' ME.message];
        errordlg(sErr);
    end
end