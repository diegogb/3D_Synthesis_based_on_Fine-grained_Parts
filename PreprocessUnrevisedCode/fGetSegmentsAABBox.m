%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2016
%Funcion: fGetSegmentsAABBox
%Input:
%   segmentsArray:  array with all the segments
%   samplingsArray: array containing the sampled 3D Point Clouds (Input
%   Shapes)
%Output:
%   segBBoxes:  array of structures containing the AABB of each segment of the array segmentsArray
%
% This function returns an array with all the Axis Aligned Bounding Boxes
% of the segments in the array "allSegments"
%
%%--------------------------------------------------------------------------
function [err, segBBoxes, segSizes, segVolumes] = fGetSegmentsAABBox(segmentsArray, samplingsArray)
    err= 0;
    segBBoxes = [];
    
    try
        numSegs= size(segmentsArray, 2);
        segSizes= zeros(numSegs, 3);
        segVolumes= zeros(numSegs, 1);
        for i = 1 : numSegs;
            bBox= fGetBoundBoxForPoints(samplingsArray(segmentsArray(i).ixShape).samples(segmentsArray(i).ixSamples, :));
            if (~isempty(bBox))
                segBBoxes= [segBBoxes bBox];
                segSizes(i, :)= bBox.extents;
                segVolumes(i)= bBox.extents(1) * bBox.extents(2) * bBox.extents(3);
            else
                err= -1;
                break;
            end
        end
        if (err)
            segBBoxes= [];
        end
    catch
        err= -1;
        segBBoxes= [];
    end
end