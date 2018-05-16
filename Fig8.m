%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to replicate the same synthesized shapes (Chairs models) presented in Figure 8 of our paper:
% 	Diego Gonzalez and Oliver van Kaick, "3D Synthesis of Man-made Objects based on Fine-grained Parts", Computers & Graphics - SMI 2018 (to appear), 2018.
%%%%%
%
%Usage: 
% 1) IMPORTANT: 
%   The variable "mainDir" below should be a valid path containing the pre-processed information 
%   of the Tables dataset. That is, "mainDir" is a folder that MUST contain the files: synthBasedOnFG-Chairs.mat and synthBasedOnFG-ChairsMS.mat
%   By default, this script assumes that the files synthBasedOnFG-Chairs.mat and synthBasedOnFG-ChairsMS.mat
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
% For example, the file "new_Fig8A-PC.ply" contains the synthesized point cloud of Figure 8(a), 
% and the file "new_Fig8A-Mesh.obj", contains the synthesized mesh of Figure 8(a)
%--------------------------------------------------------------------------

%Main Folder, which MUST contain the .mat file with the info. of the Chairs dataset:
clear;
mainDir= pwd;
% mainDir= 'C:\3D_Synthesis_based_on_Fine-grained_Parts\';

%try to find the .mat variable in the folder "mainDir"
if (exist(mainDir, 'dir'))
    if (exist(fullfile(mainDir, 'synthBasedOnFG-Chairs.mat'), 'file'))
        fprintf('This script can take a few minutes...\n');
        fprintf('Loading pre-processed data of the Chairs dataset\n');
        %Load the pre-processed information of the Chairs set:
        load(fullfile(mainDir, 'synthBasedOnFG-Chairs.mat'));
        load(fullfile(mainDir, 'synthBasedOnFG-ChairsMS.mat'));
        
        %Folder to write the generated 3D files:
        outdir= fullfile(mainDir, 'SynthesizedRep');
        if (~exist(outdir, 'dir'))
            %Create the output folder if it does not exist
            mkdir(outdir);
        end
        
        %We set the segments of the Template. For this script, all the shapes are generated from the same template
        templateIndices= fGetSegIndicesPerShape(segmentsArray, 21);
        
        fprintf('Generating shape of Figure 8(a), (and Figure 1 left)...\n');
        %%%
        %NOTE: 
        %By default, our method samples segments at random from the array of segments "segmentsArrayT", 
        %to generate new 3D shapes replacing each fine-grained segment of the template.
        %
        %For this script, we skip this sampling step, and we directly assign the segments that were 
        %selected by our method at the moment of the generation of the results shown in the paper
        selSeg= selSegsFig8a;
        %%%
        
        %Synthesize the 3D point cloud of Figure 8(a)
        [err, ~, transforms] = fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, selSeg, segmentsAdj, fullfile(outdir, 'new_Fig8A-PC.ply'), true);
        %Synthesize the 3D mesh of Figure 8(a)
        if (~err)
            err= fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSeg, transforms, templateIndices, segmentsAdj, true, fullfile(outdir, 'new_Fig8A-Mesh.obj'), false);
        end
        
        if (~err)
            %For figure 8(b)...
            fprintf('Generating shape of Figure 8(b)...\n');            
            selSeg= selSegsFig8b;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, selSeg, segmentsAdj, fullfile(outdir, 'new_Fig8B-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSeg, transforms, templateIndices, segmentsAdj, true, fullfile(outdir, 'new_Fig8B-Mesh.obj'), false);
            end
        end
        
        if (~err)
            %For figure 8(c)...
            fprintf('Generating shape of Figure 8(c)...\n');            
            selSeg= selSegsFig8c;
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, selSeg, segmentsAdj, fullfile(outdir, 'new_Fig8C-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSeg, transforms, templateIndices, segmentsAdj, true, fullfile(outdir, 'new_Fig8C-Mesh.obj'), false);
            end
        end
        
        if (~err)
            fprintf(['The process finished successfully! All the 3D Shapes of Figure 8 of the paper have been generated in the folder:\n' strrep(outdir,'\','\\')]);
            fprintf('\n');
        else
            fprintf('The script ended with errors. Some of the 3D shapes of Figure 8 of the paper probably were not generated\n');
        end
    else
        errordlg('The required files "synthBasedOnFG-Chairs.mat" and "synthBasedOnFG-ChairsMS.mat" could not be found in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the files "synthBasedOnFG-Chairs.mat" and "synthBasedOnFG-ChairsMS.mat"');
end
