%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to replicate the synthesized Table models presented in Figures 11(a) and (b) of our paper:
% 	Diego Gonzalez and Oliver van Kaick, "3D Synthesis of Man-made Objects based on Fine-grained Parts", Computers & Graphics - SMI 2018 (to appear), 2018.
%%%%%
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
% For example, the file "new_Fig11A-PC.ply" contains the synthesized point cloud of Figure 11(a), 
% and the file "new_Fig11A-Mesh.obj", contains the synthesized mesh of Figure 11(a)
%
% 3) TO SYNTHESIZE THE Chair MODELS OF FIGURE 11 OF THE PAPER, PLEASE RUN THE SCRIPT synthesizePreserveChairs_Fig11cAnd11d.M
%   FIGURES 11 (e) and (f) USE OUR HYBRID DATASET. TO SYNTHESIZE THE SHAPES FOR THESE TWO FIGURES, USE THE SCRIPT synthesizePreserveHyb_Fig11eAnd11f.M
%--------------------------------------------------------------------------

%Main Folder, which MUST contain the .mat file with the info. of the Tables dataset:
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
        
        fprintf('Generating shape of Figure 11(a)...\n');
        
        %We set the segments of the Template for the first shape: Figure 11(a)
        templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 14); 
        %%%
        %NOTE: 
        %By default, our method samples segments at random from the array of segments "segmentsArrayT", 
        %to generate new 3D shapes replacing each fine-grained segment of the template.
        %
        %To synthesize shapes preserving some segments, we keep some segments of the template, 
        %and we directly assign the rest of the segments as the ones that were selected 
        %by our method at the moment of the generation of the results shown in the paper
        selSeg= templateIndices;
        
        selSeg(7)= selSegsFig11a(7);
        selSeg(8)= selSegsFig11a(8);
        selSeg(9)= selSegsFig11a(9);
        selSeg(10)= selSegsFig11a(10);
        selSeg(12)= selSegsFig11a(12);
        selSeg(13)= selSegsFig11a(13);
        selSeg(14)= selSegsFig11a(14);
        selSeg(15)= selSegsFig11a(15);
        %%%
        
        %Synthesize the 3D point cloud of Figure 11(a)
        [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig11A-PC.ply'), true);
        %Synthesize the 3D mesh of Figure 11(a)
        if (~err)
            err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig11A-Mesh.obj'), true);
        end
        
        if (~err)
            %For figure 11(b)...
            fprintf('Generating shape of Figure 11(b)...\n');            
            
            templateIndices= fGetSegIndicesPerShape(segmentsArrayT, 8);
            selSeg= selSegsFig11b;
            selSeg(9)= templateIndices(9);
            selSeg(11)= templateIndices(11);
            selSeg(15)= templateIndices(15);
            selSeg(21)= templateIndices(21);
            selSeg(14)= templateIndices(14);
            selSeg(17)= templateIndices(17);
            selSeg(20)= templateIndices(20);
            selSeg(24)= templateIndices(24);
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsArrayT, samplingsArrayT, templateIndices, selSeg, segmentsAdjT, fullfile(outdir, 'new_Fig11B-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsT, segmentsArrayT, samplingsArrayT, selSeg, transforms, templateIndices, segmentsAdjT, false, fullfile(outdir, 'new_Fig11B-Mesh.obj'), true);
            end
        end
        
        if (~err)
            fprintf(['The process finished successfully! The 3D Shapes of Figures 11(a) and (b) of the paper have been generated in the folder:\n' strrep(outdir,'\','\\')]);
            fprintf('\n');
        else
            fprintf('The script ended with errors. Some of the 3D shapes probably were not generated\n');
        end
    else
        errordlg('The required files "synthBasedOnFG-Tables.mat" and "synthBasedOnFG-TablesMS.mat" could not be found in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the file "synthBasedOnFG-Tables.mat"');
end
