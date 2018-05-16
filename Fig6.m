%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to replicate the same synthesized Table shapes presented in Figure 6 of our paper:
% 	Diego Gonzalez and Oliver van Kaick, "3D Synthesis of Man-made Objects based on Fine-grained Parts", Computers & Graphics - SMI 2018 (to appear), 2018.
%%%%%
%
%Usage: 
% 1) IMPORTANT: 
%   The variable "mainDir" below should be a valid path containing the pre-processed information 
%   of the Tables dataset. That is, "mainDir" is a folder that MUST contain the files: synthBasedOnFG-Tables.mat and synthBasedOnFG-TablesMS.mat
%   By default, this script assumes that the files synthBasedOnFG-Tables.mat and synthBasedOnFG-TablesMS.mat
%   is in the same folder of this script.
%   %%%
%   NOTE: if desired, you can copy the ".mat" files to any directory in your computer, 
%       and then, you just have to set the variable "mainDir" to that folder
%   %%%
%
% 2) Run this script. If the process is successful, you will find one pair
% of files for each 3D shape synthesized in the folder given by the path: <mainDir>\SynthesizedRep
% One file is the synthesized point cloud: .ply file; and the other, the synthesized mesh: .obj file
%
% For example, the file "new_Fig6A-PC.ply" contains the synthesized point cloud of Figure 6(a), 
% and the file "new_Fig6A-Mesh.obj", contains the synthesized mesh of Figure 6(a)
%--------------------------------------------------------------------------

%The current folder, SHOULD contain the .mat file with the info. of the Tables dataset:
clear;
mainDir= pwd;
% mainDir= 'C:\3D_Synthesis_based_on_Fine-grained_Parts\';

%try to find the .mat variable in the folder "mainDir"
if (exist(mainDir, 'dir'))
    if (exist(fullfile(mainDir, 'synthBasedOnFG-Tables.mat'), 'file')) && (exist(fullfile(mainDir, 'synthBasedOnFG-TablesMS.mat'), 'file'))
        fprintf('This script can take a few minutes...\n');
        fprintf('Loading pre-processed data of the Tables dataset\n');
        %Load the pre-processed information of the Tables set:
        load(fullfile(mainDir, 'synthBasedOnFG-Tables.mat'));
        load(fullfile(mainDir, 'synthBasedOnFG-TablesMS.mat'));
        
        %Folder to write the generated 3D files:
        outdir= fullfile(mainDir, 'SynthesizedRep');
        if (~exist(outdir, 'dir'))
            %Create the output folder if it does not exist
            mkdir(outdir);
        end
        
        fprintf('Generating shape of Figure 6(a)...\n');
        
        %We set the segments of the Template for the first shape: Figure 6(a)
        templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 13);        
        %%%
        %NOTE: 
        %By default, our method samples segments at random from the array of segments "segmentsArrayT", 
        %to generate new 3D shapes replacing each fine-grained segment of the template.
        %
        %For this script, we skip this sampling step, and we directly assign the segments that were 
        %selected by our method at the moment of the generation of the results shown in the paper
        selSeg= selSegsFig6a;
        %%%
        
        %Synthesize the 3D point cloud of Figure 6(a)
        [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig6A-PC.ply'), true);
        %Synthesize the 3D mesh of Figure 6(a)
        if (~err)
            err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig6A-Mesh.obj'), true);
        end
        
        if (~err)
            %For figure 6(b)...
            fprintf('Generating shape of Figure 6(b)...\n');
            
            templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 10);
            selSeg= selSegsFig6b;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig6B-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig6B-Mesh.obj'), true);
            end
        end
        
        if (~err)
            %For figure 6(c)...
            fprintf('Generating shape of Figure 6(c), (and Figure 1 middle)...\n');
            templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 7);
            selSeg= selSegsFig6c;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig6C-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig6C-Mesh.obj'), true);
            end
        end
        
        if (~err)
            %For figure 6(d)...
            fprintf('Generating shape of Figure 6(d)...\n');
            templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 9);
            selSeg= selSegsFig6d;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig6D-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig6D-Mesh.obj'), true);
            end
        end
        
        if (~err)
            %For figure 6(e)...
            fprintf('Generating shape of Figure 6(e)...\n');
            templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 15);
            selSeg= selSegsFig6e;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig6E-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, true, fullfile(outdir, 'new_Fig6E-Mesh.obj'), false);
            end
        end
        
        if (~err)
            %For figure 6(f)...
            fprintf('Generating shape of Figure 6(f)...\n');
            templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 1);
            selSeg= selSegsFig6f;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig6F-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig6F-Mesh.obj'), true);
            end
        end
        
        if (~err)
            fprintf(['The process finished successfully! All the 3D Shapes of Figure 6 of the paper have been generated in the folder:\n' strrep(outdir,'\','\\')]);
            fprintf('\n');
        else
            fprintf('The script ended with errors. Some of the 3D shapes of Figure 6 of the paper probably were not generated\n');
        end
    else
        errordlg('The required files "synthBasedOnFG-Tables.mat" and "synthBasedOnFG-TablesMS.mat" could not be found in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the files "synthBasedOnFG-Tables.mat" and "synthBasedOnFG-TablesMS.mat"');
end
