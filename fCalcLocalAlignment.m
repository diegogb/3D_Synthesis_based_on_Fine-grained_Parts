%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Sep/2017 - XXX/2017
%Funcion: fCalcLocalAlignment
% Function that computes an alignment composed of a scaling and a translation, 
% to find the best possible alignment for the points "targPoints" towards
% the points "srcPoints"
% This function obtains a 3x3 scaling matrix and a translation vector, which can be transformed 
% into a homogenous transformation matrix by the calling function
%%%%%
% NOTE: this code is an adaptation of the code for alignment (translation and scaling) between two 
%   sets of points by Han Liu, which in turn, is based on code from Paper (a)
%%%%%
%(1) References: 
%   (a) Zheng, Y., et al. "Smart Variations: Functional Substructures for Part Compatibility" (2013)
%%%%%
%--------------------------------------------------------------------------
% function scalingMat = fGetScaling(segmentsAdj, shapeIndex, segIndex, contactsInfo)
function [scalingMat, transVec] = fCalcLocalAlignment(srcPoints, targPoints, maxScaling, minScCust)
    maxScFact= 1.6;
    minScFact= 1.0;
    scalingMat= [];
    try
        if (exist('maxScaling', 'var'))
            if (maxScaling < maxScFact)
                maxScFact= maxScaling;
            end
        end
        if (exist('minScCust', 'var'))
            if (minScCust > 0.5)
                minScFact= minScCust;
            end
        end
        
        numCl= size(srcPoints, 1);
        if (numCl >= 2)
            %Adapting Han's code, we compute the scaling factor per coordinate, using the info
            %of the source and targets points (vertices) of the current segment and its neighbors
            centSrc= mean(srcPoints);
            centTarg= mean(targPoints);
            
            %%%
            transVec= centTarg - centSrc;
            %%%
            
            newScX= 0;
            newScY= 0;
            newScZ= 0;
            for ixfs= 1 : numCl
                pS= srcPoints(ixfs) - centSrc;
                pT= targPoints(ixfs) - centTarg;
                newScX = newScX + (pT(1) / pS(1));
                newScY = newScY + (pT(2) / pS(2));
                newScZ = newScZ + (pT(3) / pS(3));
            end
            newScX = newScX / numCl;
            newScY = newScY / numCl;
            newScZ = newScZ / numCl;

            newScX(isnan(newScX))= 1.0;
            newScY(isnan(newScY))= 1.0;
            newScZ(isnan(newScZ))= 1.0;
            
            if (newScX <= 0)
                newScX= 1.0;
            end
            if (newScY <= 0)
                newScY= 1.0;
            end
            if (newScZ <= 0)
                newScZ= 1.0;
            end
            
%             newScX(newScX < 0.0)= minScFact;
%             newScY(newScY < 0.0)= minScFact;
%             newScZ(newScZ < 0.0)= minScFact;

            %%%
            %We try another scaling factor, and we choose the largest valid factor. 
            %This is because, in general, we're trying to enlarge the segments for filling holes.
            %NOTE: this additional factor and all the additional heuristics below, 
            %are Not from Han's original code
%             sz2= (max(targPoints)-min(targPoints)) ./ (max(srcPoints)-min(srcPoints));
% 
%             if ((sz2(1)>newScX) && (sz2(1)<maxScFact))
%                 newScX= sz2(1);
%             end
%             if ((sz2(2)>newScY) && (sz2(2)<maxScFact))
%                 newScY= sz2(2);
%             end
%             if ((sz2(3)>newScZ) && (sz2(3)<maxScFact))
%                 newScZ= sz2(3);
%             end

%             if (sz2(1) > newScX)
%                 if (sz2(1) <= maxScFact)
%                     newScX= sz2(1);
%                 else
%                     newScX= maxScFact;
%                 end
%             end
%             if (sz2(2) > newScY)
%                 if (sz2(2) <= maxScFact)
%                     newScY= sz2(2);
%                 else
%                     newScY= maxScFact;
%                 end
%             end
%             if (sz2(3) > newScZ)
%                 if (sz2(3) <= maxScFact)
%                     newScZ= sz2(3);
%                 else
%                     newScZ= maxScFact;
%                 end
%             end
             
            %%%
            %Diego Gonzalez: we add "limits" to the scaling as a way to penalize 
            %large differences w.r.t the size of the original segments
            newScX(newScX < minScFact)= minScFact;
            newScY(newScY < minScFact)= minScFact;
            newScZ(newScZ < minScFact)= minScFact;

            newScX(newScX > maxScFact)= maxScFact;
            newScY(newScY > maxScFact)= maxScFact;
            newScZ(newScZ > maxScFact)= maxScFact;
            %%%
            
            scalingMat(1,1)= newScX; 
            scalingMat(2,2)= newScY;
            scalingMat(3,3)= newScZ;
        end        
    catch
        %In case of error, we return an empty matrix 
%         scalingMat= eye(3);
        scalingMat= [];
        transVec= [];
    end    
end