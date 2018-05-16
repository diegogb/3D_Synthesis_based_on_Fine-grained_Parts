%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Feb/2016-Abr/2016
%Function: points_write_ply
%Input:
%   filename:       full name (including path) of the file where the "point cloud" will be saved
%   pointSampled:   array containing the 3D coords. of the points to write
%                   into the PLY file.
%   norms:          (optional) normals of the point sampled.
%   colors:         (optional) colors to assign to each point. (Must be an
%                   array with the RGB values)
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%
% This function saves a "point cloud" representation (points sampled from a
% Mesh), into a "PLY" file. The file includes the coords. of the points (and, optionally,
% the normals of those points)
%--------------------------------------------------------------------------
% IMPORTANT NOTE: This function is based on the script "mesh_write_ply.m" developed by, 
% Oliver van Kaick, Professor (and Diego Gonzalez's Thesis Advisor) at Carleton University.
%--------------------------------------------------------------------------
function [err] = points_write_ply(filename, pointSampled, norms, colors)
    err = 0;
    incNorms = 0;
    incCols= 0;
    fid = -1;
    try
        % Open file in write mode
        fid = fopen(filename, 'w');
        if (fid ~= -1)
            fprintf(fid, 'ply\n');
            fprintf(fid, 'format ascii 1.0\n');
            fprintf(fid, 'element vertex %d\n', size(pointSampled, 1));
            fprintf(fid, 'property float x\n');
            fprintf(fid, 'property float y\n');
            fprintf(fid, 'property float z\n');
            
            if exist('colors', 'var') && ~isempty(colors)
                incCols = 1;
                fprintf(fid, 'property uchar red\n');
                fprintf(fid, 'property uchar green\n');
                fprintf(fid, 'property uchar blue\n');
            end
            
            if exist('norms', 'var') && ~isempty(norms)
                incNorms = 1;
                fprintf(fid, 'property float nx\n');
                fprintf(fid, 'property float ny\n');
                fprintf(fid, 'property float nz\n');
            end
            fprintf(fid, 'end_header\n');

%             % Write points info into the file
%             if (incCols && incNorms)
%                 fprintf(fid, '%f %f %f %d %d %d %f %f %f\n', [pointSampled(:, 1) pointSampled(:, 2) pointSampled(:, 3) ...
%                     uint8(colors(:, 1)) uint8(colors(:, 2)) uint8(colors(:, 3)) norms(:, 1) norms(:, 2) norms(:, 3)]);
%             else
%                 if (incCols)
%                     fprintf(fid, '%f %f %f %d %d %d\n', [pointSampled(:, 1) pointSampled(:, 2) pointSampled(:, 3) ...
%                         uint8(colors(:, 1)) uint8(colors(:, 2)) uint8(colors(:, 3))]);
%                 elseif(incNorms)
%                     fprintf(fid, '%f %f %f %f %f %f\n', [pointSampled(:, 1) pointSampled(:, 2) pointSampled(:, 3) ...
%                         norms(:, 1) norms(:, 2) norms(:, 3)]);
%                 else
%                     fprintf(fid, '%f %f %f\n', [pointSampled(:, 1) pointSampled(:, 2) pointSampled(:, 3)]);
%                 end
%             end
            
            for i = 1:size(pointSampled, 1)
%             for i = 1:size(pointSampled.point, 1)
                %1st: Point coords.
                fprintf(fid, '%f %f %f', pointSampled(i, 1), pointSampled(i, 2), pointSampled(i, 3));
%                 fprintf(fid, '%f %f %f', pointSampled.point(i, 1), pointSampled.point(i, 2), pointSampled.point(i, 3));
                if (incCols)
                     %2nd: colors of the points.
                    fprintf(fid, ' %d %d %d', uint8(colors(i, 1)), uint8(colors(i, 2)), uint8(colors(i, 3)));
                end
                if (incNorms)
                     %3rd: normals of the points.
                    fprintf(fid, ' %f %f %f', norms(i, 1), norms(i, 2), norms(i, 3));
                end
                fprintf(fid, '\n');
            end

            % Close file
            fclose(fid);
        else
            sErr= 'Error saving a PLY file. The file could not be open for writing';
            errordlg(sErr);
        end    
    catch ME
        if (fid ~= -1)
            % Close file
            fclose(fid);
        end
        sErr= ['Error saving a PLY file (error in points_write_ply): ' ME.message];
        errordlg(sErr);
    end
end