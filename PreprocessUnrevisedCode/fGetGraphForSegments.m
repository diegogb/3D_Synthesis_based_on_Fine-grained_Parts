%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2016
%Funcion: fGetGraphForSegments
%Input:
%   allSegments:    structure array containing the segments that were previously computed for all 
%                   the input shapes
%   samplingsArray: array of structures containing the sampled points from the Meshes
%   ixShape:        index of the Shape for which the graph will be computed
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%   segAdjacency:   pairwise adjacency Matrix indicating the segments of
%                   the shape "shapeName" that are neighbors
%   adjInfo:        structure containing additional info. about the
%                   neighbouring segments: main axis indicated by the min.
%                   variance of the contact points
%
% This function obtains a graph (an adjacency Matrix) for all the segments of the Shape "ixShape"
%
% (See also fGetAdjacencyPerShape, fGetSegIndicesPerShape)
%
%%%%%
%(1) References: 
%   (a) Prakhar, J., et al. "Assembly-based conceptual 3D modeling with unlabeled components using probabilistic factor graph" (2016)
%%%%%
%--------------------------------------------------------------------------
function [err, segAdjacency, adjInfo] = fGetGraphForSegments(allSegments, samplingsArray, ixShape)
% function [err, segAdjacency, adjInfo] = fGetGraphForSegments(allSegments, samplingsArray, ixShape, segSymArray)
    err= 0;
    
    volThres= 0.035;
    try
        %1st. find the indices of the segments of the Shape "shapeName"
        if (size(samplingsArray, 2) > 0) && (ixShape <= size(samplingsArray, 2))
            segIndices= fGetSegIndicesPerShape(allSegments, ixShape);
        
            if (~isempty(segIndices))
                numSegsShape= length(segIndices);
                contacts= struct([]);
                adjInfo= struct([]);
                defMainAx= zeros(numSegsShape, 1);
%                 segMainAx= zeros(numSegsShape, 1);
                tmpMainAx= zeros(numSegsShape, 3);

                %Get the volume of the Shape:
                szShape= max(samplingsArray(ixShape).samples) - min(samplingsArray(ixShape).samples);
                volShape = szShape(1)*szShape(2)*szShape(3);

%                 thres= 1/32;
                
%                 %We use a threshold based on a percentage of the Shape's
%                 %volume, as in Paper (a)
                thres= (volShape^(1/3)) * volThres;
%                 thres= (volShape*(1/3)) * volThres;
%                 thres= 1/32;
                
                segAdjacency= zeros(numSegsShape, numSegsShape);

                %We check every pair of segments of the Shape "ixShape" to determine the neighbour segments
                for i= 1 : numSegsShape
%                     %We first use PCA as an initial way to learn the global axis 
%                     %that is most similar to the most important orientation of the segment
                    ptsI= samplingsArray(ixShape).samples(allSegments(segIndices(i)).ixSamples, :);
                    szSegI= max(ptsI) - min(ptsI);
                    [~, tempA]= max(szSegI);
                    centI= mean(ptsI);
                    
%                     mC= cov(ptsI);
%                     [eVecs, ~]= eigs(mC);
% 
%                     vS= eVecs(:, 3)';
%                     tempA= 0;
%                     maxDp= -1;
%                     for a= 1 : 3
%                         aX= zeros(1, 3);
%                         aX(a)= 1;
%                         dp= abs(sum(vS .* aX));
%                         if (dp > maxDp)
%                             maxDp= dp;
%                             tempA= a;
%                         end
%                     end

                    %We use the axis of larger extent as initial default orientation: 
                    %this is only a heuristic that we try to improve below
                    %using information about the neighbouring segments
                    defMainAx(i)= tempA;
                    tmpMainAx(i, tempA)= tmpMainAx(i, tempA) + 1;

                    for j= i + 1 : numSegsShape
                        contacts(i, j).ixSample= [];
                        contacts(i, j).coord= [];
                        contacts(j, i).ixSample= [];
                        contacts(j, i).coord= [];

                        adjInfo(i, j).parAx= [];
                        adjInfo(i, j).side= [];
                        %Get the distances between each pair of points of each
                        %segment
                        ptsJ= samplingsArray(ixShape).samples(allSegments(segIndices(j)).ixSamples, :);
                        distBS= pdist2(ptsI, ptsJ);
                        
%                         numClose= length(unique(find(distBS(:) <= thres)));
%                         if (numClose >= 4)
                        isClose= find(distBS(:) <= thres, 1);
                        if (~isempty(isClose))
                            %We consider two segments i, j adjacent, if there's at least  a pair of points, 
                            %one in each segment, closer than our given threshold. (Reference (a))
                            
                            segAdjacency(i, j)= 1;
                            segAdjacency(j, i)= 1;

                            %We extract 10 closest points of each segment and
                            %we store them as "contact" points                        
                            [~, ixSor]= sort(distBS(:));
                            [closeI, closeJ]= ind2sub(size(distBS), ixSor);
                            numCI= 1;
                            ixClose= 1;
                            while (ixClose<= 10 && numCI< size(closeI, 1))
                                %Store the contacts for the "i" segment
                                if (~ismember(closeI(numCI), contacts(i, j).ixSample))
                                    contacts(i, j).ixSample= cat(1, contacts(i, j).ixSample, closeI(numCI));
                                    contacts(i, j).coord= cat(1, contacts(i, j).coord, ...
                                        samplingsArray(ixShape).samples(allSegments(segIndices(i)).ixSamples(closeI(numCI)), :));

                                    ixClose= ixClose + 1;                                
                                end
                                numCI= numCI + 1;
                            end
%                             centContI= mean(contacts(i, j).coord);
                            
                            numCJ= 1;
                            ixClose= 1;
                            while (ixClose<= 10 && numCJ< size(closeJ, 1))
                                %Store the contacts for the "j" segment
                                if (~ismember(closeJ(numCJ), contacts(j, i).ixSample))
                                    contacts(j, i).ixSample= cat(1, contacts(j, i).ixSample, closeJ(numCJ));
                                    contacts(j, i).coord= cat(1, contacts(j, i).coord, ...
                                        samplingsArray(ixShape).samples(allSegments(segIndices(j)).ixSamples(closeJ(numCJ)), :));

                                    ixClose= ixClose + 1;
                                end                            
                                numCJ= numCJ + 1;                            
                            end
%                             centContJ= mean(contacts(j, i).coord);

                            %%%
                            %We store additional info about the adjacent segment: global axis that is 
                            %most similar to the vector between the neighbours
                            %centroids, and sign (0,-,+) of the subtraction of those centroids                            
                            centJ= mean(samplingsArray(ixShape).samples(allSegments(segIndices(j)).ixSamples, :));

                            [~, tempA]= min(var(contacts(i, j).coord));

                            adjInfo(i, j).parAx= tempA;
                            adjInfo(i, j).side= sign(centI(tempA) - centJ(tempA));

                            tmpMainAx(i, tempA)= tmpMainAx(i, tempA) + 1;

                            [~, tempA]= min(var(contacts(j, i).coord));

                            adjInfo(j, i).parAx= tempA;
                            adjInfo(j, i).side= -adjInfo(i, j).side;

                            tmpMainAx(j, tempA)= tmpMainAx(j, tempA) + 1;
                            
%                             vBS= centI - centJ;
%                             vBS= vBS / norm(vBS);
%                             tempA= 0;
%                             maxDp= -1;
%                             for a= 1 : 3
%                                 aX= zeros(1, 3);
%                                 aX(a)= 1;
%                                 dp= abs(sum(vBS .* aX));
%                                 if (dp > maxDp)
%                                     maxDp= dp;
%                                     tempA= a;
%                                 end
%                             end
%                             
%                             tmpMainAx(i, tempA)= tmpMainAx(i, tempA) + 1;
%                             tmpMainAx(j, tempA)= tmpMainAx(j, tempA) + 1;
                            %%%
                        end
                    end
                end
                
%                 %TODO: Check the graph's symmetry: symmetric segments
%                 %should have a consistent number of segments, (not always
%                 %the same number, but it should be consistent)
%                 segChk= false(numSegsShape, 1);
%                 for i= 1 : numSegsShape
%                     if (~segChk(i))
%                         segChk(i)= true;
%                         %Get the symmetric segment: we get only the "first one"
%                         candiDel= [];
%                         prevC= [];
%                         [err, symSeg]= fGetSymmetryChain(segIndices(i), allSegments, segIndices);                    
%     %                     symSeg= allSegments(segSymArray(segIndices(i)).ixSymSegments(iSy)).segLabel;
%                         if (~err && ~isempty(symSeg))
%                             segChk(symSeg)= true;
%                             symSeg= symSeg(symSeg ~= i);
%                             for iSy= 1 : length(symSeg)
%                                 ixNeighs= find(segAdjacency(i, :) == 1);
%                                 ixNeighsSymm= find(segAdjacency(symSeg(iSy), :) == 1);
% 
%                                 %We'll keep the smallest number of neighbors found between the two symmetric segs.
%                                 if (length(ixNeighs) > length(ixNeighsSymm)) || (~isempty(candiDel))
%                                     for nh= 1 : length(ixNeighs)
%                                         [~, symGNeigh] = fGetSymmetryChain(segIndices(ixNeighs(nh)), allSegments, segIndices);
%                                         if (~isempty(symGNeigh))
%     %                                         symGNeigh= symGNeigh(symGNeigh ~= ixNeighs(nh));
%                                             for iSyS= 1 : length(symGNeigh)
%             %                                     symGNeigh= allSegments(segSymArray(segIndices(ixNeighs(nh))).ixSymSegments).segLabel;
%                                                 if (~ismember(symGNeigh(iSyS), ixNeighsSymm)) && (~ismember(ixNeighs(nh), prevC))
%                                                     candiDel= cat(1, candiDel, ixNeighs(nh));
%             %                                                 segAdjacency(i, ixNeighs(nh))= 0;
%             %                                                 segAdjacency(ixNeighs(nh), i)= 0;
%                                                 else
%                                                     if (ismember(ixNeighs(nh), candiDel))
%                                                         prevC= cat(1, prevC, ixNeighs(nh));
%                                                         candiDel= candiDel(candiDel ~= ixNeighs(nh));
%                                                     end
%                                                     break;
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                         candiDel= unique(candiDel);
%                         for d= 1 : length(candiDel)
%                             segAdjacency(i, candiDel(d))= 0;
%                             segAdjacency(candiDel(d), i)= 0;
%                         end
%                     end
%                 end
            else
                err= -1;
                errordlg('Error obtaining the graph for the segments of a Shape. Check if the index of the Template Shape is valid');
            end
        else
            err= -1;
            errordlg('Error obtaining the graph for the segments of a Shape, (error in the input parameters for "fGetGraphForSegments")');
        end        
    catch ME
        err= -1;        
        errordlg(['Error obtaining the graph for the fine-grained segments of a Shape (error in "fGetGraphForSegments"): ' ME.message]);
    end
end