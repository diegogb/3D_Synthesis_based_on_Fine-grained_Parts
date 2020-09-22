%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: May/2016
%Function: fGetFarthestPointsFromSeg
%Input:
%   points:         set of 3D points that will be sampled
%   minNumSamples:  (OPTIONAL) minimum No. of samples for "small" segments.
%                   If not given, this number is computed as: max(100, 0.1*NumPoints)
%Output:
%   err:        -1 if some error ocurrs; 0, otherwise
%   sampledIxs: array containing the indices of the points sampled
% 
% This function performs a Farthest Point (sub)Sampling of the points of a
% segment, (given by the input parameter "points"). The number of samples 
% can be given by the parameter "minNumSamples", or it will be set by
% default to "max(100, 0.1*NumPoints)"
%
% (See also: fCalcSDFHistoForSegment)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Kaick, Oliver Van, Fish, N., et al. "Shape segmentation by approximate convexity analysis." ACM (TOG) 34.1. (2014)
%   (b) Mémoli, F., and Sapiro, G. "Comparing point clouds". Procs. of Eurographics/ACM SIGGRAPH Symposium on Geom. Proc. ACM, (2004).
%(2) This code is based on Reference (b) and the implementation developed by NOA FISH for Paper (a).
%   Code available online at: http://www.cs.tau.ac.il/~noafish/wcseg/
%(3) Do Not use the parameter "minNumSamples"  to use the default number of samples: max(100, 0.1*NumPoints)
%(4) If the input has less than 100 points, all the points will be returned, 
%   assuming that the subsampling is not required (useful) in this case
%%%%%
%--------------------------------------------------------------------------
function [err, sampledIxs] = fGetFarthestPointsFromSeg(points, minNumSamples)
    err= 0;
    try
        if (size(points, 2) == 3)
            numPts= size(points, 1);
            if (~exist('minNumSamples', 'var'))
                %Set the number of samples, following Paper (a)
                minNumSamples= floor(max(100, 0.1*numPts));
            end
            
            if (numPts > minNumSamples)
                %Get the distances between every pair of points in the segment
                allD= pdist2(points, points);

                doSel = false(1, numPts);
                min_distances = inf(1, numPts);
                %We always start with the first point of the segment by default
                selection = 1;
                for i = 1 : minNumSamples
                    %Save the selected point
                    doSel(selection) = true;

                    selection_distances = zeros(1, numPts);
                    notSel = ~doSel;
                    %Get the distance between the currently selected point and all
                    %the rest of points that have not been selected
                    selection_distances(notSel)= allD(selection, notSel);
                    
                    %Sort the distances, to get the farthest point w.r.t
                    %the points that are already selected
                    min_distances = min(min_distances, selection_distances);
                    [~,srtd] = sort(min_distances, 'descend');
                    selection = srtd(1);
                end
                %We could return the logical matrix of selected points, but we return
                %the indices of the points instead
                sampledIxs= find(doSel);
            else
                sampledIxs= 1 : size(points, 1);            
            end
        else
            err= -1;
            errordlg('Error in "fGetFarthestPointsFromSeg". The input parameters are not valid');
        end
    catch ME
        err= -1;
        errordlg(['Error sampling farthest points from a segment (error in "fGetFarthestPointsFromSeg"): ' ME.message]);
    end
end