%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2017-XXX/2018
%Funcion: fSynthesizeNewShapeCollectionP
%Input:
%   segmentsArray:  structure array containing the segments that were previously computed for all 
%                   the input shapes
%   samplingsArray: array of structures containing the sampled points from the Meshes
%   meshesArraySD:  array containing the subdivided input meshes (see fGetSubDividedMeshes)
%   Dist:           Matrix of distances or similarity between each pair of
%                   segments (see fGetDistanceBetweenSegs)
%   segmentsAdj:    structure array with one element per shape, containing a graph containing the
%                   connectivity between segments, and information about
%                   the main axis of orientation of each segment (see fGetAdjacencyPerShape)
%   scoresArray:    array of scores computed for each segment, according to
%                   Eq. (5) in our paper, (see fCalcScoresWithNeighbours)
%   statsPerShape:  array with the mean and std dev. of each shape, 
%                   computed from the values of the "scoresArray". (see fCalcScoresWithNeighbours)
%   minProbThres:	minimum probability threshold for selecting the
%                   replacement of each segment of the Template. Required value is between 0.1 and 1.0, (0.75 by default)
%   savingPath:     folder to write the files (.obj and .ply) for the synthesized shapes
%Output:
%   err:                -1 if some error ocurrs; 0, otherwise
%   selSegC:            matrix containing the indices of the segments selected
%                       as replacements of the original segments of the template shape
%   selSegScoreC:       matrix with the Energies of each replacement segment,
%                       computed
%   newSegmentsArray:   (OPTIONAL) complete array of the new fine-grained segments used
%                       to synthesize each 3D shape
%                       NOTE: this can be redundant information
%   avgTimeSynthPC:     (OPTIONAL) time required by the process of Synthesizing a Point cloud, 
%                       including writing a .ply file with the point cloud info.
%   avgTimeSynthMesh:   (OPTIONAL) time required by the process of Synthesizing a Mesh, 
%                       including writing an .obj file with the mesh info.
% Function that generates a new shape collection, using each one of shapes stored
% in the structure array "samplingsArray" as a Template.
%
% (See also fGetAdjacencyPerShape, fGetSubDividedMeshes, fGetDistanceBetweenSegs, 
%    fCalcScoresWithNeighbours, fSampleSegsForNewShape, fSynthesizeNovelShape, fSynthesizeMeshFromPoints)
%
%--------------------------------------------------------------------------
function [err, selSegC, selSegScoreC, newSegmentsArray, avgTimeSynthPC, avgTimeSynthMesh] = ...
    fSynthesizeNewShapeCollectionP(segmentsArray, samplingsArray, meshSegments, Dist, segmentsAdj, scoresArray, statsPerShape, minProbThres, savingPath)
    err= 0;    
    try
        if (~exist('minProbThres', 'var'))
            minProbThres= 0.75;
        end
        
        if (minProbThres < 0.1 || minProbThres > 1.0)
            0.75;
        end
        
        %Synthesize a novel shape using as template each Shape of the input
        %dataset
        numS= size(samplingsArray, 2);        
        selSegC= zeros(numS, 80);
        selSegScoreC= zeros(numS, 80);
        newSegmentsArray= struct([]);
        %%%
        %To store processing times:
        timePC= zeros(1, numS);
        timeMesh= zeros(1, numS);
        %%%
        
        %For each shape in the Input dataset...
        for ixt= 1 : numS
            %Get the indices of the segments of the current template shape:
            templateIndices= fGetSegIndicesPerShape(segmentsArray, ixt);
            
            fprintf(['Sampling new segments to synthesize a 3D shape with Template No.' num2str(ixt,'%d') '...\n']);
            %Sample segments that will replace the segments of the template
            [err, selSeg, selSegScore] = fSampleSegsForNewShape(segmentsArray, samplingsArray, segmentsAdj, templateIndices, Dist, ...
                scoresArray, statsPerShape, minProbThres, false);
            if (err)
                break;
            end

            selSegC(ixt, 1:length(selSeg))= selSeg;
            selSegScoreC(ixt, 1:length(selSeg))= selSegScore;
            thePath= fullfile([savingPath 'new_' num2str(ixt,'%d') '_' num2str(sum(selSegScore), '%.4f')]);
            
            %Now, synthesize the point cloud for the given template, using
            %the sampled segments
            fprintf(['Synthesizing point cloud with Template No.' num2str(ixt,'%d') '...\n']);
            tic;
            [err, newSegments, transforms]= fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, ...
                selSeg, segmentsAdj, fullfile([thePath '.ply']));
            timePC(ixt)= toc;
            if (~err)
                fprintf(['Synthesizing mesh with Template No.' num2str(ixt,'%d') '...\n']);
                tic;
                err= fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSeg, transforms, templateIndices, segmentsAdj, ...
                    true, fullfile([thePath '.obj']), false);
                timeMesh(ixt)= toc;
            else
                break;
            end
            
            newSegmentsArray= [newSegmentsArray newSegments];
        end
        
        avgTimeSynthPC= mean(timePC);        
        avgTimeSynthMesh= mean(timeMesh);
    catch ME
        err= -1;        
        errordlg(['Error synthesizing a collection of shapes (error in "fSynthesizeNewShapeCollection"): ' ME.message]);
    end
end