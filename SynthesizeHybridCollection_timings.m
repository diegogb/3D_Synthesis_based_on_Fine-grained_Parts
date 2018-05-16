%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%%%%%
% This script allows you to generate a new collection of 3D shapes from our input set of Tables, as well as a text file with timing statistics for
% the synthesis process using our approach, which is presented in:
% 	Diego Gonzalez and Oliver van Kaick, "3D Synthesis of Man-made Objects based on Fine-grained Parts", Computers & Graphics - SMI 2018 (to appear), 2018.
%%%%%
%
%Usage:
%Usage: 
% 1) IMPORTANT: 
%   The variable "mainDir" below should be a valid path containing the pre-processed information 
%   of the Hybrid dataset. That is, "mainDir" is a folder that MUST contain the 3 files: 
%       synthBasedOnFG-Hybrid.mat, synthBasedOnFG-HybridMS1.mat, synthBasedOnFG-HybridMS2.mat
%   By default, this script assumes that the files
%   synthBasedOnFG-Hybrid.mat, synthBasedOnFG-HybridMS1.mat, and synthBasedOnFG-HybridMS2.mat
%   are in the same folder of this script.
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
%
%
% %%%
% NOTES ABOUT THE TIMING:
% %%%
% Additionally, you will find in the folder <mainDir>\Synthesized a file called "timings.txt" with two times:
% 1) Average Time to Synthesize a Point cloud
% 2) Average Time to Synthesize a Mesh
% ******    ******
%   IMPORTANT: 
%   1) If the file "timings.txt" already exists, this script will append the new times at the end of the file
%   2) The times reported in our article, were obtained as the average of four executions of this script
% ******    ******
% %%%
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
    if (exist(fullfile(mainDir, 'synthBasedOnFG-Hybrid.mat'), 'file')) && (exist(fullfile(mainDir, 'synthBasedOnFG-HybridMS1.mat'), 'file'))
        fprintf('This script can take several minutes...\n');
        fprintf('Loading pre-processed data of the Hybrid dataset\n');
        %Load the pre-processed information of the Hybrid set:
        load(fullfile(mainDir, 'synthBasedOnFG-Hybrid.mat'));
        load(fullfile(mainDir, 'synthBasedOnFG-HybridMS1.mat'));
        load(fullfile(mainDir, 'synthBasedOnFG-HybridMS2.mat'));
        
        %We had to divide the data of the segments of the hybrid in two
        %parts, so here we concat. the info. again
        meshSegmentsAll= [meshSegmentsA1 meshSegmentsA2];
        
        %Folder to write the generated 3D files:
        outdir= fullfile(mainDir, 'Synthesized');
        if (~exist(outdir, 'dir'))
            %Create the output folder if it does not exist
            mkdir(outdir);
        end
        
        [err, ~, ~, ~, avgTimeSynthPC, avgTimeSynthMesh] = ...
            fSynthesizeNewShapeCollectionP(segmentsAll, samplingsAll, meshSegmentsAll, DistAll, segmentsAdjAll, ...
            scoresArrayAll, statsPerShapeAll, taoProbThres, [outdir '\']);
        
        if (~err)
            fid = -1;
            try
                % Open timings file in append mode
                filename= fullfile(outdir, 'timings.txt');
                fid = fopen(filename, 'a');
                if (fid ~= -1)
                    curT= clock;
                    curTstr= [num2str(curT(1)) '-' num2str(curT(2)) '-' num2str(curT(3)) ' ' num2str(curT(4)) ':' num2str(curT(5)) ':' num2str(curT(6), 2)];
                    msgToP= ['Hybrid set: synthesis of shape collection ending at:' curTstr];
                    msgToP= [msgToP ' - Average time (seconds) to synthesize each point cloud: %.4f\n'];
                    fprintf(fid, msgToP, avgTimeSynthPC);
                    
                    msgToP= ['Hybrid set: synthesis of shape collection ending at:' curTstr];
                    msgToP= [msgToP ' - Average time (seconds) to synthesize each mesh: %.4f\n'];
                    fprintf(fid, msgToP, avgTimeSynthMesh);
                    
                    % Close file
                    fclose(fid);
                    
                    fprintf(['The process finished successfully! A new collection of 3D models created with the hybrid set is in the folder:\n' strrep(outdir,'\','\\')]);
                    fprintf('\n');
                    fprintf('This folder also contains the file "timings.txt" with the average time used to synthesize each shape');
                    fprintf('\n');
                else
                    sErr= 'Error saving the file "timings.txt". The file could not be opened for writing';
                    errordlg(sErr);
                end    
            catch ME
                if (fid ~= -1)
                    % Close file
                    fclose(fid);
                end    

                errordlg(ME.message);
            end            
        end
    else
        errordlg('The required ".mat" files of the hybrid dataset could not be found in the folder given by the variable "mainDir"');
    end
else
    errordlg('The variable "mainDir" must specify a valid folder containing the required .mat files of the hybrid set');
end
