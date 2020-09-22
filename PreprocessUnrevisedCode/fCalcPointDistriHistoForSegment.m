%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Jul/2017 - XXX/2017
%Function: fCalcPointDistriHistoForSegment
%Input:
%   points:         set of 3D points for which the Histogram of the Point
%                   distribution will be computed, following Paper (a)
%Output:
%   err:        -1 if some error ocurrs; 0, otherwise
%   pdistHisto: histogram of size 64 for the set of points given as input
% 
% This function computes a histogram that estimates the distribution of the
% points, by overlaying a grid of a give size (4 x 4 x 4), and counting the
% number of points inside each cell of the grid
% 
%
% (See also fGetDistanceBetweenSegs, fGetSegmentsAABBox, fGetBoundBoxForPoints)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Lun, Z., et al. "Functionality preserving shape style transfer." ACM (TOG) 35.6 (2016).
%%%%%
%--------------------------------------------------------------------------
function [err, pdistHisto] = fCalcPointDistriHistoForSegment(points)
    %Initialize return value with zeros
    pdistHisto= zeros(64, 1);
    err= 0;
    numGridCells= 4;
    try
        if (size(points, 1) > 1 && size(points, 2) == 3)
            numP= size(points, 1);
            
            %Get the bounding box for the points
            bBox= fGetBoundBoxForPoints(points);
            
            %Get the size of the grid per axis
            gridSize= bBox.extents ./ numGridCells;

            %Get grids
            gX= min(points(:, 1)) : gridSize(1) : max(points(:, 1));
            gY= min(points(:, 2)) : gridSize(2) : max(points(:, 2));
            gZ= min(points(:, 3)) : gridSize(3) : max(points(:, 3));
            
            for ix= 1 : numGridCells
                for iy= 1 : numGridCells
                    for iz= 1 : numGridCells
                        if (ix < numGridCells)
                            insideX= find(points(:, 1) >= gX(ix) & points(:, 1) < gX(ix+1));
                        else
                            insideX= find(points(:, 1) >= gX(ix) & points(:, 1) <= gX(ix+1));
                        end
                        if (iy < numGridCells)
                            insideY= find(points(:, 2) >= gY(iy) & points(:, 2) < gY(iy+1));
                        else
                            insideY= find(points(:, 2) >= gY(iy) & points(:, 2) <= gY(iy+1));
                        end
                        if (iz < numGridCells)
                            insideZ= find(points(:, 3) >= gZ(iz) & points(:, 3) < gZ(iz+1));
                        else
                            insideZ= find(points(:, 3) >= gZ(iz) & points(:, 3) <= gZ(iz+1));
                        end
                        
                        inBin= intersect(intersect(insideX, insideY), insideZ);
                        if (~isempty(inBin))
                            pdistHisto(ix*iy*iz)= size(inBin, 1);
                        end
                    end
                end
            end
            
            %Normalize dividing over the number of points
            if (~err)
                pdistHisto= pdistHisto / numP;
            end            
        else
            err= -1;
            errordlg('Error in fCalcPointDistriHistoForSegment: invalid input sizes');
        end
    catch ME
        err= -1;
        errordlg(['Error computing a histogram for the point distribution of a segment (error in fCalcPointDistriHistoForSegment): ' ME.message]);
    end
end