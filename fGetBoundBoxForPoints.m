%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2016
%Funcion: fGetBoundBoxForPoints
%Input:
%   points:    matrix of size [N*3] containing the point set for which
%              the Bounding Box will be computed
%Output:
%   bBox:   structure containing each corner or vertex of the Bounding
%           box, and the centroid of the points
%
% This function returns the bounding box of a set of 3D points, given by
% the vertices (corners) of the box, its centroid (mean point), and the
% extents (size) of the box
% (If there is an error, the function returns an "empty" value= []
%
%%--------------------------------------------------------------------------
function bBox = fGetBoundBoxForPoints(points)
    bBox= [];
    try
        if (size(points, 2) == 3)
            bBox.center = mean(points);
            tempMax= max(points);
            
            bBox.v1= min(points);
            bBox.v2= [bBox.v1(1) bBox.v1(2) tempMax(3)];
            bBox.v3= [bBox.v1(1) tempMax(2) bBox.v1(3)];
            bBox.v4= [bBox.v1(1) tempMax(2) tempMax(3)];     
            bBox.v5= [tempMax(1) bBox.v1(2) bBox.v1(3)];
            bBox.v6= [tempMax(1) bBox.v1(2) tempMax(3)];
            bBox.v7= [tempMax(1) tempMax(2) bBox.v1(3)];
            bBox.v8= tempMax;
            
            bBox.extents= tempMax - bBox.v1;
            
            %Faces of the Bounding Box:
            facets= zeros(6, 4, 3);
            facets(1,:,:) = cat(1, bBox.v1, bBox.v2, bBox.v3, bBox.v4);
            facets(2,:,:) = cat(1, bBox.v5, bBox.v6, bBox.v7, bBox.v8);
            facets(3,:,:) = cat(1, bBox.v3, bBox.v4, bBox.v7, bBox.v8);
            facets(4,:,:) = cat(1, bBox.v1, bBox.v2, bBox.v5, bBox.v6);
            facets(5,:,:) = cat(1, bBox.v1, bBox.v3, bBox.v5, bBox.v7);
            facets(6,:,:) = cat(1, bBox.v2, bBox.v4, bBox.v6, bBox.v8);

            bBox.faces= facets;
        else
            bBox= [];
            errordlg('Error: the function "fGetBoundBoxForPoints" requires as input an "N * 3" matrix (a 3D point set)');
        end        
    catch ME
        bBox= [];
        errordlg(['Error computing the Bounding Box (error in fGetBoundBoxForPoints): ' ME.message]);
    end
end