function generate_pbs_mat(pbsFolderBaseName, pbsFolderCoreName, videoFileName, frames, templateFileName, outFileName, coilSensitivityFile, hpcFolder, funFolder)
%GENERATE_PBS Summary of this function goes here
%   Detailed explanation goes here

pbsFolderName = sprintf('%s/%s',pbsFolderBaseName,pbsFolderCoreName);

mkdir(pbsFolderName);
mkdir([pbsFolderName,'/experiment/']);
mkdir([pbsFolderName,'/logs/']);

copyfile([funFolder,'/calckbkernel.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/switch_contour_tracker.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/gridkb.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/gridlut.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/ift.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/interpft2.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/kb.m'],[pbsFolderName,'/experiment']);
copyfile([funFolder,'/gridlut_mex.mexa64'],[pbsFolderName,'/experiment']);


copyfile(templateFileName,[pbsFolderName, '/experiment/template.mat'])

if ~isempty(coilSensitivityFile)
   copyfile(coilSensitivityFile,[pbsFolderName,'/experiment/coilSensitivity.mat']); 
end;

videostruct=avi_to_struct(videoFileName); videodata=videostruct.frames;
save([pbsFolderName,'/experiment/videodata.mat'],'videodata');

fid=fopen([pbsFolderName,'/experiment/main.m'],'w');
fprintf(fid, 'distcomp.feature( ''LocalUseMpiexec'', true );\n');
fprintf(fid, 'load template.mat;\n');
fprintf(fid, 'load videodata.mat;\n');
fprintf(fid, 'frames = %s;\n',mat2str(frames));
fprintf(fid, 'numIterations = [10 90 300 300];\n');
if ~isempty(coilSensitivityFile)
    fprintf(fid, 'coilIntensityCorrectionOn = 1;\n');
else
    fprintf(fid, 'coilIntensityCorrectionOn = 0;\n');
end
fprintf(fid, 'newtonMethodOn = 0;\n');
fprintf(fid, 'contourCleanUpOn = 1;\n');
fprintf(fid, 'plotOn = 0;\n');
fprintf(fid, 'notifyOn = 0;\n');
fprintf(fid, 'parallelOn = 1;\n');
fprintf(fid, 'workers = 8;\n');
fprintf(fid, 'coilSensitivityMatFileName = ''coilSensitivity.mat'';\n');

fprintf(fid, 'trackdata = switch_contour_tracker(videodata, template_struct, numIterations,...\n');
fprintf(fid, '    coilIntensityCorrectionOn, coilSensitivityMatFileName, ...\n');
fprintf(fid, '    newtonMethodOn, contourCleanUpOn, plotOn, notifyOn, frames, ...\n');
fprintf(fid, '    parallelOn, workers);\n');

fprintf(fid, 'save(''../../%s_track.mat'',''trackdata'')\n',outFileName);

fclose(fid);

fid=fopen(sprintf('%s/pbs',pbsFolderName),'w');

fprintf(fid, '#PBS -l walltime=23:59:59\n');
%fprintf(fid, '#PBS -l mem=1gb\n');
fprintf(fid, '#PBS -l nodes=1:ppn=8\n');
fprintf(fid, '#PBS -o %s/%s/logs/output_main.txt\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, '#PBS -e %s/%s/logs/error_main.txt\n\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, 'source /usr/usc/matlab/default/setup.sh\n');
fprintf(fid, 'cd %s/%s/experiment\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, 'export MATLAB_PREFDIR=%s/%s/scratch\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, 'source /usr/usc/matlab/default/setup.sh\n');
%fprintf(fid, 'sleep $[ ( $RANDOM %% 120 )  + 1 ]s\n');
fprintf(fid, 'matlab -r \"main;exit;\"\n');

fclose(fid);

end

