%--------------------------------------------------------------------------
%%%%%
% (Diego Gonzalez): This function is MODIFIED for changing the required parameter
% to be a set of 3D points (representing a point cloud), instead of a Mesh.
% The paramter "samples" contains points sampled from a Shape, (a part of a 3D Shape)
%Input:
%   samples:        set of 3D points for which the OBB will be computed
%   sortToGlobal:   (OPTIONAL) Indicates if the axes will be sorted
%                   according to the global Shape's axes, as described in
%                   Reference (a)
%Output:
%   err:        -1 if some error ocurrs; 0, otherwise
%%%%%
% NOTES:
% This function was originally written by Noa Fish, and it is modified here for working
% with a set of points (a point cloud), instead of requiring a triangular mesh
% THE ORIGINAL VERSION OF THIS CODE WAS WRITTEN FOR THE Paper:
%   (a) Fish, N., Averkiou, M., van Kaick, O., et al. "Meta-representation of shape families". (2014)
%%%%%
%--------------------------------------------------------------------------
function bbox = part_obb(samples, sortToGlobal)
% function bbox = part_obb(samples, mesh)
    bbox= [];
    try
        thresh = 1e-4;
        
        if (~exist('sortToGlobal', 'var'))
            sortToGlobal= true;
        end
        
        %%%
%         faces = mesh.faces;
%         points = mesh.vertices;
%         coh = convhull(points(:,1), points(:,2), points(:,3));
%         faces= [];
%         points = samples;
        coh = convhull(samples(:,1), samples(:,2), samples(:,3));
        %%%
        
        currDir= pwd();
        newDir= fullfile(currDir, 'obb');
        if (exist(newDir, 'dir'))
            cd(newDir);
            
            %%%
            % (Diego Gonzalez):
            % tic;
            %[mbox, boxes] = point_sym_oobb(samples, points, coh, faces);
%             [mbox, boxes] = sym_bb(samples, points, coh, faces);
            % The function symOBB_PointSet is a C++ function that makes all the actual work of computing the OBB
            [mbox, boxes] = symOBB_PointSet(samples, coh, 0);
            %%%
            % toc;
            %%%
            [sbox, axes, ~] = find_sym_box(boxes, thresh);

            %sym = [];
            rot = [];

            if(isempty(sbox))
                box = mbox;
                bx = process_bbox(box);
            else
                box = reshape(sbox(1:end-3),3,5);
                bx = process_bbox(box);

                %sym = axes;

                %%%
                % (Diego Gonzalez): rotational symmetry is not being checked
                % check rotational sym
                %rot = check_sym2(box, axes, mesh, samples, thresh);
                %%%
            end

            bbox = struct('center', bx(:,1), 'origin', bx(:,2), 'axes', [box(:,2:4); 2*box(:,5)'], 'sym', axes, 'rot', rot);
            cd(currDir);
            
            bbox = calc_axis_dir(bbox);
        end
        
        %%%
        %Diego Gonzalez: we add the sorting of the OBB axes to best align
        %them to the global axes, following Reference (a)
        if (~isempty(bbox) && sortToGlobal)
            ordAx= zeros(1, 3);
            dp= abs(sum(bbox.axes(1:3, :) .* repmat([1 0 0]', 1, 3)));
%             dp= sum(bbox.axes(1:3, :) .* repmat([1 0 0]', 1, 3));
            [~, ordAx(1)]= max(dp);
            dp= abs(sum(bbox.axes(1:3, :) .* repmat([0 1 0]', 1, 3)));
            [~, ordAx(2)]= max(dp);
            dp= abs(sum(bbox.axes(1:3, :) .* repmat([0 0 1]', 1, 3)));
            [~, ordAx(3)]= max(dp);
            
            %Check if there is more than 1 axis of the OBB pairing with each global axis
            if (length(unique(ordAx, 'stable')) == 3)
                bbox.axes= bbox.axes(:, ordAx);
            else                
                dp= sum(bbox.axes(1:3, :) .* repmat([1 0 0]', 1, 3));
                [~, ordM]= max(dp);
                if (ordM ~= ordAx(1))
                    ordAx(1)= ordM;
                end
                dp= sum(bbox.axes(1:3, :) .* repmat([0 1 0]', 1, 3));
                [~, ordM]= max(dp);
                if (ordM ~= ordAx(2))
                    ordAx(2)= ordM;
                end
                dp= sum(bbox.axes(1:3, :) .* repmat([0 0 1]', 1, 3));
                [~, ordM]= max(dp);
                if (ordM ~= ordAx(3))
                    ordAx(3)= ordM;
                end
                
                if (length(unique(ordAx)) == 3)
                    bbox.axes= bbox.axes(:, ordAx);
                end
            end
        end

        %Diego Gonzalez: also, we add to the object that is returned a matrix
        %representing the faces of the OBB
        %Get the corners of the OBB of the reference segment
        center = bbox.center;
        exts = bbox.axes(end, :);
        hexts = 0.5*exts;
        axs = bbox.axes(1:3,:);
        origin = center - axs * hexts';

        v11 = bbox.axes(1:3, 1) * bbox.axes(4, 1);
        v22 = bbox.axes(1:3, 2) * bbox.axes(4, 2);
        v33 = bbox.axes(1:3, 3) * bbox.axes(4, 3);    
        p1 = origin;
        p2 = p1 + v11;
        p3 = p1 + v22;
        p4 = p1 + v33;
        p5 = p2 + v22;
        p6 = p2 + v33;
        p7 = p6 + v22;
        p8 = p7 - v11;

        %Faces of the OBB of the current segment:
        facets= zeros(6, 4, 3);
        facets(1,:,:) = [p1 p2 p5 p3]';
        facets(2,:,:) = [p1 p2 p6 p4]';
        facets(3,:,:) = [p1 p3 p8 p4]';
        facets(4,:,:) = [p2 p5 p7 p6]';
        facets(5,:,:) = [p3 p5 p7 p8]';
        facets(6,:,:) = [p4 p6 p7 p8]';
        
        bbox.faces= facets;
        %%%
    catch ME
        bbox= [];
        sErr= ['Error computing the OBB of a fine-grained segment (error in part_obb): ' ME.message];
        errordlg(sErr);        
    end
    
end
%%%%%   %%%%%   %%%%%
% END Changes By Diego Gonzalez
%%%%%   %%%%%   %%%%%

function [sbox, axes, bid] = find_sym_box(boxes, thresh)
    %thresh = 1e-5; %1e-3;

    sbox = [];
    axes = [];
    bid = -1;
    %rots = [];

    num = size(boxes,2);
    min1 = 10;
    min2 = 10;
    min3 = 10;
    min1_ind = -1;
    min2_ind = -1;
    min3_ind = -1;

    for i=1:num
        ext1 = boxes(end-5,i);
        ext2 = boxes(end-4,i);
        ext3 = boxes(end-3,i);
        vol = ext1 * ext2 * ext3;
        score1 = boxes(end-2,i);
        score2 = boxes(end-1,i);
        score3 = boxes(end,i);

        if(vol < min1 && score1 < thresh)
            min1 = vol;
            min1_ind = [i 1];
        end

        if(vol < min1 && score2 < thresh)
            min1 = vol;
            min1_ind = [i 2];
        end

        if(vol < min1 && score3 < thresh)
            min1 = vol;
            min1_ind = [i 3];
        end

        if(vol < min2 && score1 < thresh && score2 < thresh)
            %min2 = (score1+score2)/2;
            min2 = vol;
            min2_ind = [i 1 2];
        end

        if(vol < min2 && score1 < thresh && score3 < thresh)
            %min2 = (score1+score3)/2;
            min2 = vol;
            min2_ind = [i 1 3];
        end

        if(vol < min2 && score2 < thresh && score3 < thresh)
            %min2 = (score2+score3)/2;
            min2 = vol;
            min2_ind = [i 2 3];
        end

        if(vol < min3 && score1 < thresh && score2 < thresh && score3 < thresh)
            %min3 = (score1+score2+score3)/3;
            min3 = vol;
            min3_ind = i;
        end

    end

    %return;

    if(min3_ind ~= -1)
         sbox = boxes(:,min3_ind);
         axes = 1:3;
         bid = min3_ind;
    %      rot1 = is_rot_sym(1, sbox, boxes, thresh);
    %      rot2 = is_rot_sym(2, sbox, boxes, thresh);
    %      rot3 = is_rot_sym(3, sbox, boxes, thresh);
    %      
    %      rots = [rot1, rot2, rot3];
    %      nz = find(rots > 0);
    %      rots = rots(nz);

         return;
    end

    if(min2_ind ~= -1)
         sbox = boxes(:,min2_ind(1));
         axes = min2_ind(2:3);
         bid = min2_ind(1);
    %      rot_cand = setdiff(1:3, axes);
    %      rot1 = is_rot_sym(rot_cand, sbox, boxes, thresh);
    %      
    %      rots = rot1;
    %      nz = find(rots > 0);
    %      rots = rots(nz);

         return;
    end

    if(min1_ind ~= -1)
         sbox = boxes(:,min1_ind(1));
         axes = min1_ind(2);
         bid = min1_ind(1);
         return;
    end
end

function bx = process_bbox(box)
    cen = box(:,1);
    ext = box(:,5);
    [sr_ext, sr_ind] = sort(ext, 'descend');
    v1 = box(:,1+sr_ind(1));
    v2 = box(:,1+sr_ind(2));
    v3 = box(:,1+sr_ind(3));


    v11 = v1*2*sr_ext(1);
    v22 = v2*2*sr_ext(2);
    v33 = v3*2*sr_ext(3);

    %origin = cen - v1*ext(1) - v2*ext(2) - v3*ext(3);
    origin = cen - v1*sr_ext(1) - v2*sr_ext(2) - v3*sr_ext(3);

    bx = [cen origin v11 v22 v33];

end


