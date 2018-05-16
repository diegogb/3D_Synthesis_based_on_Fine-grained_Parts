%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to replicate the synthesized Chair models presented in Figures 11(c) and (d) of our paper:
% 	Diego Gonzalez and Oliver van Kaick, "3D Synthesis of Man-made Objects based on Fine-grained Parts", Computers & Graphics - SMI 2018 (to appear), 2018.
%%%%%
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
% For example, the file "new_Fig11C-PC.ply" contains the synthesized point cloud of Figure 11(c), 
% and the file "new_Fig11C-Mesh.obj", contains the synthesized mesh of Figure 11(c)
%
% 3) TO SYNTHESIZE THE Tables MODELS OF FIGURE 11 OF THE PAPER, PLEASE RUN THE SCRIPT synthesizePreserveTables_Fig11aAnd11b.M
%   FIGURES 11 (e) and (f) USE OUR HYBRID DATASET. TO SYNTHESIZE THE SHAPES FOR THESE TWO FIGURES, USE THE SCRIPT synthesizePreserveHyb_Fig11eAnd11f.M
%--------------------------------------------------------------------------

%Main Folder, which MUST contain the .mat file with the info. of the Chairs dataset:
clear;
mainDir= pwd;
% mainDir= 'C:\3D_Synthesis_based_on_Fine-grained_Parts\';

%try to find the .mat variable in the folder "mainDir"
if (exist(mainDir, 'dir'))
    if (exist(fullfile(mainDir, 'synthBasedOnFG-Chairs.mat'), 'file'))
        fprintf('This script can take a few minutes...\n');
        fprintf('Loading pre-processed data of the Chair dataset\n');
        %Load the pre-processed information of the Chairs set:
        load(fullfile(mainDir, 'synthBasedOnFG-Chairs.mat'));
        load(fullfile(mainDir, 'synthBasedOnFG-ChairsMS.mat'));
        
        %Folder to write the generated 3D files:
        outdir= fullfile(mainDir, 'SynthesizedRep');
        if (~exist(outdir, 'dir'))
            %Create the output folder if it does not exist
            mkdir(outdir);
        end
        
        fprintf('Generating shape of Figure 11(c)...\n');
        
        %We set the segments of the Template for the first shape: Figure 11(c)
        templateIndices= fGetSegIndicesPerShape(segmentsArray, 4); 
        %%%
        %NOTE: 
        %By default, our method samples segments at random from the array of segments "segmentsArray", 
        %to generate new 3D shapes replacing each fine-grained segment of the template.
        %
        %To synthesize shapes preserving some segments, we keep some segments of the template, 
        %and we directly assign the rest of the segments as the ones that were selected 
        %by our method at the moment of the generation of the results shown in the paper
        selSeg= selSegsFig11c;
        
        selSeg(1)= templateIndices(1);
        selSeg(2)= templateIndices(2);
        selSeg(3)= templateIndices(3);
        selSeg(4)= templateIndices(4);
        %%%
        
        %Synthesize the 3D point cloud of Figure 11(c)
        [err, ~, transforms] = fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, selSeg, segmentsAdj, fullfile(outdir, 'new_Fig11C-PC.ply'), true);
        %Synthesize the 3D mesh of Figure 11(c)
        if (~err)
            err= fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSeg, transforms, templateIndices, segmentsAdj, true, fullfile(outdir, 'new_Fig11C-Mesh.obj'), false);
        end
        
        if (~err)
            %For figure 11(b)...
            fprintf('Generating shape of Figure 11(d)...\n');
            
            templateIndices= fGetSegIndicesPerShape(segmentsArray, 10);
            selSeg= selSegsFig11d;
            selSeg(45)= templateIndices(45);
            selSeg(46)= templateIndices(46);
            selSeg(47)= templateIndices(47);
            selSeg(48)= templateIndices(48);
            selSeg(56)= templateIndices(56);
            
            %Synthesize the 3D point cloud of Figure 11(d)
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArray, samplingsArray, templateIndices, selSeg, segmentsAdj, fullfile(outdir, 'new_Fig11D-PC.ply'), true);
            %Synthesize the 3D mesh of Figure 11(d)
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegments, segmentsArray, samplingsArray, selSeg, transforms, templateIndices, segmentsAdj, false, fullfile(outdir, 'new_Fig11D-Mesh.obj'), true);
            end
        end
        
        if (~err)
            fprintf(['The process finished successfully! The 3D Shapes of Figures 11(c) and (d) of the paper have been generated in the folder:\n' strrep(outdir,'\','\\')]);
            fprintf('\n');
        else
            fprintf('The script ended with errors. Some of the 3D shapes probably were not generated\n');
        end
    else
        errordlg('The variable "synthBasedOnFG-Chairs.mat" could not be found. Please check that this variable is in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the file "synthBasedOnFG-Chairs.mat"');
end
