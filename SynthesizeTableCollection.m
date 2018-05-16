%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to generate a new collection of 3D shapes from our input set of Tables 
% using our approach. See details in our paper:
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
% 2) We recomend keeping the default Probability threshold "tao" of 0.75
% (see our paper for details). You can, however, perform experiments
% changing this value: variable "taoProbThres"
% NOTE: the valid range for "taoProbThres" is between 0.1 and 1. However,
% using values taoProbThres<= 0.5 could generate highly implausible 3D shapes
%
% 3) Run this script. If the process is successful, you will find one pair
% of files for each 3D shape synthesized in the folder given by the path: <mainDir>\Synthesized
% One file is the synthesized point cloud: .ply file; and the other, the synthesized mesh: .obj file
% ******
%   IMPORTANT NOTE ABOUT FILE NAMES: 
%       The files names are generated automatically, using the prefix "new"
%       and two additional "parts":
%       1) a unique number for each Template Shape of the Tables dataset
%       2) The Shape Energy computed with Eq. (2) of our paper.
%       EXAMPLES:
%       1) File "new_1_95.8933.ply": Point cloud synthesized for the template No. 1, which has a Shape Energy= "95.8933"
%       2) File "new_1_95.8933.obj": 3D mesh synthesized for the template No. 1, which has a Shape Energy= "95.8933"
% ******
%--------------------------------------------------------------------------
%The current folder, SHOULD contain the .mat file with the info. of the Tables dataset:
clear;
mainDir= pwd;
% mainDir= 'C:\3D_Synthesis_based_on_Fine-grained_Parts\';

%Probability threshold tao to sample segments from the Probability
%Distribution previously generated (see comments above and our paper for details)
taoProbThres= 0.75;

%try to find the .mat variable in the folder "mainDir"
if (exist(mainDir, 'dir'))
    if (exist(fullfile(mainDir, 'synthBasedOnFG-Tables.mat'), 'file')) && (exist(fullfile(mainDir, 'synthBasedOnFG-TablesMS.mat'), 'file'))
        fprintf('This script can take several minutes...\n');
        fprintf('Loading pre-processed data of the Tables dataset\n');
        %Load the pre-processed information of the Tables set:
        load(fullfile(mainDir, 'synthBasedOnFG-Tables.mat'));
        load(fullfile(mainDir, 'synthBasedOnFG-TablesMS.mat'));
        
        %Folder to write the generated 3D files:
        outdir= fullfile(mainDir, 'Synthesized');
        if (~exist(outdir, 'dir'))
            %Create the output folder if it does not exist
            mkdir(outdir);
        end
        
        [err, ~, ~, ~, avgTimeSynthPC, avgTimeSynthMesh] = ...
            fSynthesizeNewShapeCollectionP(segmentsArrayT, samplingsArrayT, meshSegmentsT, DistT, segmentsAdjT, ...
            scoresArrayT, statsPerShapeT, taoProbThres, [outdir '\']);
        
        if (~err)
            fprintf(['The process finished successfully! A new collection of 3D table models is in the folder:\n' strrep(outdir,'\','\\')]);
            fprintf('\n');
        end
    else
        errordlg('The required files "synthBasedOnFG-Tables.mat" and "synthBasedOnFG-TablesMS.mat" could not be found in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the files "synthBasedOnFG-Tables.mat" and "synthBasedOnFG-TablesMS.mat"');
end
