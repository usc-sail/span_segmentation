addpath(genpath(fullfile(pwd,'functions')));

pbsBaseFolderName='../cluster';
% This is the folder on your local computer where the files to upload to the cluster will be put

hpcFolder='/home/rcf-proj2/tjs/tsorense';
% This is the folder on the cluster where you will upload the data (make
% sure you have access to that folder)

mkdir(pbsBaseFolderName);
segmentFile = '../manual_annotations/timestamps.csv';
% Your .csv file
firstLine = 2;
lastLine = 13;
% The first and last lines of the part of the .csv file you need to work
% for this experiment

templateFileName = '../template_struct_converted.mat';
coilSensitivityFile = [];

frameRate = 83.33;
% You can find this as videostruct.framerate from previous steps of the
% process

scriptName = 'do_segmentation';
% This is just the name of the shell script that will be created

funFolder = fullfile(pwd,'functions/contour_tracker');
% Location when the contour_tracker functions are stored

%generate_pbs_from_file_full_mat(pbsBaseFolderName, scriptName, segmentFile, firstLine, lastLine, frameRate, templateFileName, coilSensitivityFile, hpcFolder, funFolder)
cluster_generate_scripts(pbsBaseFolderName, scriptName, segmentFile, firstLine, lastLine, frameRate, templateFileName, coilSensitivityFile, hpcFolder, funFolder)

% OK. Now go to <pbsBaseFolderName>, upload everything from there to
% hpc-transfer.usc.edu
% Then log to hpc-login3.usc.edu and execute the <scriptName>.sh script
% Hopefully, after a few hours you will get a bunch of .mat files at your
% home folder on hpc-login3.usc.edu which you can transfer back to your
% computer.