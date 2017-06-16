function generate_pbs(pbsFolderBaseName, pbsFolderCoreName, videoFileName, frames, templateFileName, outFileName, coilSensitivityFile, hpcFolder)
%GENERATE_PBS Summary of this function goes here
%   Detailed explanation goes here

pbsFolderName = sprintf('%s/%s',pbsFolderBaseName,pbsFolderCoreName);

mkdir(pbsFolderName);
mkdir([pbsFolderName,'/experiment/']);
mkdir([pbsFolderName,'/logs/']);

copyfile('./functions/contour_tracker/calckbkernel.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/switch_contour_tracker.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/gridkb.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/gridlut.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/ift.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/interpft2.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/kb.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/avi_to_struct.m',[pbsFolderName,'/experiment']);
copyfile('./functions/contour_tracker/gridlut_mex.mexa64',[pbsFolderName,'/experiment']);


copyfile(templateFileName,[pbsFolderName, '/experiment/template.mat'])
%copyfile([videodataFileName,'.mat'],[pbsFolderName,'/experiment/videodata.mat']);
copyfile(videoFileName,[pbsFolderName,'/experiment/videofile.avi']);

if ~isempty(coilSensitivityFile)
   copyfile(coilSensitivityFile,[pbsFolderName,'/experiment/coilSensitivity.mat']); 
end;

fid=fopen([pbsFolderName,'/experiment/main.m'],'w');
fprintf(fid, 'load template.mat;\n');
fprintf(fid, 'videostruct=avi_to_struct(''videofile.avi''); videodata=videostruct.frames;\n');
%fprintf(fid, 'load videodata.mat;\n');
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
fprintf(fid, 'workers = 16;\n');
fprintf(fid, 'coilSensitivityMatFileName = ''coilSensitivity.mat'';\n');

fprintf(fid, 'trackdata = switch_contour_tracker(videodata, template_struct, numIterations,...\n');
fprintf(fid, '    coilIntensityCorrectionOn, coilSensitivityMatFileName, ...\n');
fprintf(fid, '    newtonMethodOn, contourCleanUpOn, plotOn, notifyOn, frames, ...\n');
fprintf(fid, '    parallelOn, workers)\n');

fprintf(fid, 'save(''../../%s_track.mat'',''trackdata'')\n',outFileName);

fclose(fid);

fid=fopen(sprintf('%s/pbs',pbsFolderName),'w');

fprintf(fid, '#PBS -l walltime=23:59:59\n');
%fprintf(fid, '#PBS -l mem=1gb\n');
fprintf(fid, '#PBS -l nodes=1:ppn=16\n');
fprintf(fid, '#PBS -o %s/%s/logs/output_main.txt\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, '#PBS -e %s/%s/logs/error_main.txt\n\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, 'source /usr/usc/matlab/default/setup.sh\n');
fprintf(fid, 'cd %s/%s/experiment\n',hpcFolder,pbsFolderCoreName);
fprintf(fid, 'matlab -r \"main;exit;\"\n');

fclose(fid);

end

