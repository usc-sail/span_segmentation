addpath(genpath('../../../a_functions'));

pbsBaseFolderName='../cluster';
% This is the folder on your local computer where the files to upload to the cluster will be put

hpcFolder='/home/rcf-proj2/sn/toutios';
% This is the folder on the cluster where you will upload the data (make
% sure you have access to that folder)

mkdir(pbsBaseFolderName);
segmentFile = './csvfile_eff.csv';
% Your .csv file
firstLine = 1;
lastLine = 6;
% The first and last lines of the part of the .csv file you need to work
% for this experiment

templateFileName = './template_struct_converted.mat';
coilSensitivityFile = []; 

scriptName = 'split_script_1_6';
% This is just the name of the shell script that will be created

funFolder = '/Users/toutios/Dropbox/Programming/a_functions/contour_tracker';
% Location when the contour_tracker functions are stored

cluster_generate_scripts_split(pbsBaseFolderName, scriptName, segmentFile, firstLine, lastLine, templateFileName, coilSensitivityFile, hpcFolder, funFolder);

% OK. Now go to <pbsBaseFolderName>, upload everything from there to
% hpc-transfer.usc.edu
% Then log to hpc-login3.usc.edu and execute the <scriptName>.sh script
% Hopefully, after a few hours you will get a bunch of .mat files at your
% home folder on hpc-login3.usc.edu which you can transfer back to your
% computer.