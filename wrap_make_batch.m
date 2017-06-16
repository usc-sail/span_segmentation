addpath(genpath('functions'))

pbsBaseFolderName=fullfile(cd,'cluster');
% This is the folder on your local computer where the files to upload to the cluster will be put

hpcFolder='/home/rcf-proj2/tjs/tanner';
% This is the folder on the cluster where you will upload the data (make
% sure you have access to that folder)

mkdir(pbsBaseFolderName);
segmentFile = 'demo_files/segments_ms.csv';
% Your .csv file
% For converting from frames to ms:
% tab = readtable('../csv/at1_rep.csv','Delimiter',','), tab.Var3 = tab.Var3*1000/frameRate; tab.Var4 = tab.Var4*1000/frameRate; writetable(tab,'../csv/at1_rep.csv')

firstLine = 2;
lastLine = 50;
% The first and last lines of the part of the .csv file you need to work
% for this experiment

templateFileName = 'template_struct_converted.mat';
% This is the output of wrap_template_batch.m

coilSensitivityFile = [];

frameRate = 83.2778;
% You can find this as videostruct.framerate from previous steps of the
% process

scriptName = 'demo';
% This is just the name of the shell script that will be created

account = 'lc_tjs';
% This is the account whose core hours will be used.

generate_pbs_from_file(pbsBaseFolderName, scriptName, segmentFile, firstLine, lastLine, frameRate, templateFileName, coilSensitivityFile, hpcFolder, account)

% OK. Now go to <pbsBaseFolderName>, upload everything from there to
% hpc-transfer.usc.edu
% Then log to hpc-login3.usc.edu and execute the <scriptName>.sh script
% Hopefully, after a few hours you will get a bunch of .mat files at your
% home folder on hpc-login3.usc.edu which you can transfer back to your
% computer.