%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2016 - XXX/2017
%Function: fCalcAnglesForPFH
%Input:
%   points:     set of 3D points for which the Point Feauture Histogram
%               (PFH) will be computed, following Papers (a) and (b)
%   normals:    matrix of size [N * 3] with the normals of the 3D points
%Output:
%   err:      -1 if some error ocurrs; 0, otherwise
%   pfhHisto: histogram of size 125, i.e. Point Feature Histogram (PFH) for the
%             set of points and the normals given as inputs
% 
% This function computes a histogram that estimates the difference between the normals
% of each pair of points, following References (a) and (b)
%
% (See also fGetDistanceBetweenSegs)
%
%%%%%
% NOTES: 
%(1) References: 
%   (a) Zhao, X., et al. "Indexing 3D Scenes Using the Interaction Bisector Surface". (2014)
%   (b) Rusu, R., et al. "Learning informative point classes for the acquisition of object model maps". (2008)
%%%%%
%--------------------------------------------------------------------------
function [err, pfhHisto] = fCalcAnglesForPFH(points, normals)
    %Initialize return value with zeros
    pfhHisto= zeros(125, 1);
    err= 0;
    zeroEpsil= 1e-9;
    try
        if (size(points, 1) == size(normals, 1)) && (size(points, 2) == 3 && size(normals, 2) == 3)
            numP= size(points, 1);
            
            %Compute the distance between each pair of points
            dist = pdist2(points, points);

            %The angles are computed here between each point "i" and all the rest of the points of the segment
            for i= 1 : numP - 1
                %u = normal vector of current point
                u= normals(i, :);
                %Compute the difference between each pair of points in the set: (Pj - Pi)
                diffP= bsxfun(@minus, points(i+1:end, :), points(i, :));
                %Next, compute: v= u Cross ((Pi - Pj) / d) ; where u = normal of Pi
                v= cross(repmat(u, numP - i, 1), bsxfun(@rdivide, diffP, dist(i+1:end, i)));
                %Then, compute: w= u Cross v
                w= cross(repmat(u, numP - i, 1), v);

                %Store n2= normals of "2nd" point: Pj
                n2= normals(i+1:end, :);

                %Compute the 3 angles (alpha, phi, theta), following Paper (a)
                %alpha= dot(v, n2)
                alpha= sum(v .* n2, 2);
                alpha(alpha >= -zeroEpsil & alpha <= zeroEpsil)= 0;
                alpha(isnan(alpha))= 0;

                %phi= dot(u, ((Pj - Pi) / Distance_between_Points))
                phi= sum(repmat(u, numP - i, 1) .* bsxfun(@rdivide, diffP, dist(i+1:end, i)), 2);
                phi(phi >= -zeroEpsil & phi <= zeroEpsil)= 0;
                phi(isnan(phi))= 0;

                %theta= arctan(dot(w, n2), dot(u, n2))
                theta= atan2(sum(w .* n2, 2), sum(repmat(u, numP - i, 1) .* n2, 2));
    %             theta= atan(sum(w .* n2, 2) ./ sum(repmat(u, numP - i, 1) .* n2, 2));
                theta(theta >= -zeroEpsil & theta <= zeroEpsil)= 0;
                theta(isnan(theta))= 0;

                %%%
                %Now, compute the Histogram of length 125 for the 3 angles:            
                %BASED ON Reference (b), get the index of the bin where the point pair should be
                %Here, we compute first the index for each angle, in 5 bins, and later, we get the "total" index for the 125-bins histogram

                %1) For alpha and phi the range of values is [-1..1] = 2, so we divide by 2
                ix1= ceil(((alpha + 1) * 5) * 0.5);
                ix1(ix1 > 5)= 5;
                ix1(ix1 < 1)= 1;
                ix2= ceil(((phi + 1) * 5) * 0.5);
                ix2(ix2 > 5)= 5;
                ix2(ix2 < 1)= 1;
                %2) For theta, the range is [-pi .. pi]= 2 * pi, so we divide by (2 * pi)
                ix3= ceil(((theta + pi) * 5) / (2 * pi));
                ix3(ix3 > 5)= 5;
                ix3(ix3 < 1)= 1;

                %Now, we "combine" the indices of each angle
                idxH = ix1 + (ix2 * 5) + (ix3 * 25);            
                idxH(idxH > 125)= 125;
                idxH(idxH < 1)= 1;

                uniBins= unique(idxH);
                for ixb= 1 : length(uniBins)
                    pfhHisto(uniBins(ixb))= pfhHisto(uniBins(ixb)) + sum(idxH == uniBins(ixb));
                end
                %%%

    %             if any(idxH < 1 | idxH > 125)
    %                 err= -1;
    %                 pfhHisto= zeros(125, 1);
    %                 break;
    %             else
    %                 uniBins= unique(idxH);
    %                 for ixb= 1 : length(uniBins)
    %                     pfhHisto(uniBins(ixb))= pfhHisto(uniBins(ixb)) + sum(idxH == uniBins(ixb));
    %                 end
    %             end
                %%%
            end

            %Normalize dividing over the number of point pairs, (based on Paper (b))
            if (~err)
                pfhHisto= pfhHisto / (numP * (numP + 1) / 2);
            end
            
        else
            err= -1;
            errordlg('Error in fCalcAnglesForPFH: the size of the set of input Points and the size of input Normals should be the same');
        end
    catch ME
        err= -1;
        errordlg(['Error computing the Point Feature Histogram for a segment (error in fCalcAnglesForPFH): ' ME.message]);
    end
end