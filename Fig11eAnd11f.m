%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to replicate the synthesized Hybrid models presented in Figures 11(e) and (f) of our paper:
% 	Diego Gonzalez and Oliver van Kaick, "3D Synthesis of Man-made Objects based on Fine-grained Parts", Computers & Graphics - SMI 2018 (to appear), 2018.
%%%%%
%
%Usage: 
% 1) IMPORTANT: 
%   The variable "mainDir" below should be a valid path containing the pre-processed information 
%   of the Hybrid dataset. That is, "mainDir" is a folder that MUST contain the file: synthBasedOnFG-Hybrid.mat
%   By default, this script assumes that the file synthBasedOnFG-Hybrid.mat
%   is in the same folder of this script.
%   %%%
%   NOTE: if desired, you can copy the file synthBasedOnFG-Hybrid.mat to any directory in your computer, 
%       and you just have to set the variable "mainDir" to that folder
%   %%%
%
% 2) Run this script. If the process is successful, you will find one pair
% of files for each 3D shape synthesized in the folder given by the path: <mainDir>\SynthesizedRep
% One file is the synthesized point cloud: .ply file; and the other, the synthesized mesh: .obj file
%
% For example, the file "new_Fig11E-PC.ply" contains the synthesized point cloud of Figure 11(e), 
% and the file "new_Fig11E-Mesh.obj", contains the synthesized mesh of Figure 11(e)
%
% 3) TO SYNTHESIZE THE Table MODELS OF FIGURE 11 OF THE PAPER, PLEASE RUN THE SCRIPT synthesizePreserveTables_Fig11aAnd11b.m
%   TO SYNTHESIZE THE Chair MODELS OF FIGURE 11 OF THE PAPER, USE THE SCRIPT synthesizePreserveChairs_Fig11cAnd11d.m
%--------------------------------------------------------------------------

%Main Folder, which MUST contain the .mat file with the info. of the Hybrid dataset:
clear;
mainDir= pwd;
% mainDir= 'C:\3D_Synthesis_based_on_Fine-grained_Parts\';

%try to find the .mat variable in the folder "mainDir"
if (exist(mainDir, 'dir'))
    if (exist(fullfile(mainDir, 'synthBasedOnFG-Hybrid.mat'), 'file'))
        fprintf('This script can take a few moments...\n');
        fprintf('Loading pre-processed data of the Hybrid dataset\n');
        %Load the pre-processed information of the Hybrid set:
        load(fullfile(mainDir, 'synthBasedOnFG-Hybrid.mat'));
        
        %Folder to write the generated 3D files:
        outdir= fullfile(mainDir, 'SynthesizedRep');
        if (~exist(outdir, 'dir'))
            %Create the output folder if it does not exist
            mkdir(outdir);
        end
        
        fprintf('Generating shape of Figure 11(e)...\n');
        
        %We set the segments of the Template for the first shape: Figure 11(a)
        templateIndices= fGetSegIndicesPerShape(segmentsAll, 2); 
        %%%
        %NOTE: 
        %By default, our method samples segments at random from the array of segments "segmentsAll", 
        %to generate new 3D shapes replacing each fine-grained segment of the template.
        %
        %To synthesize shapes preserving some segments, we keep some segments of the template, 
        %and we directly assign the rest of the segments as the ones that were selected 
        %by our method at the moment of the generation of the results shown in the paper
        selSeg= selSegsFig11e;
        
        selSeg(1)= templateIndices(1);
        selSeg(2)= templateIndices(2);
        selSeg(3)= templateIndices(3);
        selSeg(4)= templateIndices(4);
        %%%
        
        %Synthesize the 3D point cloud of Figure 11(e)
        [err, ~, transforms] = fSynthesizeNovelShape(segmentsAll, samplingsAll, templateIndices, selSeg, segmentsAdjAll, fullfile(outdir, 'new_Fig11E-PC.ply'), true);
        %Synthesize the 3D mesh of Figure 11(e)
        if (~err)
            err= fSynthesizeMeshFromPoints(meshSegmentsAll, segmentsAll, samplingsAll, selSeg, transforms, templateIndices, segmentsAdjAll, false, fullfile(outdir, 'new_Fig11E-Mesh.obj'), true);
        end
        
        if (~err)
            %For figure 11(f)...
            fprintf('Generating shape of Figure 11(f)...\n');            
            
            templateIndices= fGetSegIndicesPerShape(segmentsAll, 18);
            selSeg= selSegsFig11f;
            selSeg(3)= templateIndices(3);
            selSeg(55)= templateIndices(55);
            selSeg(56)= templateIndices(56);
            selSeg(61)= templateIndices(61);
            selSeg(62)= templateIndices(62);
            
            [err, ~, transforms] = fSynthesizeNovelShape(segmentsAll, samplingsAll, templateIndices, selSeg, segmentsAdjAll, fullfile(outdir, 'new_Fig11F-PC.ply'), true);
            if (~err)
                err= fSynthesizeMeshFromPoints(meshSegmentsAll, segmentsAll, samplingsAll, selSeg, transforms, templateIndices, segmentsAdjAll, true, fullfile(outdir, 'new_Fig11F-Mesh.obj'), false);
            end
        end
        
        if (~err)
            fprintf(['The process finished successfully! The 3D Shapes of Figures 11(e) and (f) of the paper have been generated in the folder:\n' strrep(outdir,'\','\\')]);
            fprintf('\n');
        else
            fprintf('The script ended with errors. Some of the 3D shapes probably were not generated\n');
        end
    else
        errordlg('The variable "synthBasedOnFG-Hybrid.mat" could not be found. Please check that this variable is in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the file "synthBasedOnFG-Hybrid.mat"');
end
