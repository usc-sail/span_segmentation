function cluster_generate_scripts_update(localFolderName, scriptName, segmentFile, firstLine, lastLine, frameRate, templateFileName, coilSensitivityFile, hpc_folder, funFolder)
lineCount = 0;
frames =[];
scriptFileName = sprintf('%s/%s.m', localFolderName, scriptName);
%scriptFileName2 = sprintf('%s/%s.sh', localFolderName, scriptName);

fidin = fopen(segmentFile);
fidout = fopen(scriptFileName,'w');
%fprintf(fidout,'MDCSprofile = parcluster(''TwentyFourHr_MultiNodeProfile'');\n');


while lineCount<firstLine - 1
    
    fgetl(fidin);
    lineCount = lineCount + 1;
    
end;

while lineCount<lastLine
    
    str = fgetl(fidin);
    celldata = textscan(str,'%s','Delimiter',',');
    strdata = celldata{1};
    videoFileName = char(strdata(1));
    pbsFolderCoreName = char(strdata(2));
    outFileName = pbsFolderCoreName;
    
    pbsFolderName = sprintf('%s/%s',localFolderName,pbsFolderCoreName);
    
    mkdir(pbsFolderName);
    mkdir([pbsFolderName,'/experiment/']);
    
    copyfile([funFolder,'/calckbkernel.m'],[pbsFolderName,'/experiment']);
    copyfile([funFolder,'/cluster_contour_tracker_update.m'],[pbsFolderName,'/experiment']);
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
    fprintf(fid, 'warning(''off'',''all'');\n');
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
    %fprintf(fid, 'workers = 63;\n');
    fprintf(fid, 'coilSensitivityMatFileName = ''coilSensitivity.mat'';\n');
    
    fprintf(fid, 'try\n');
    fprintf(fid, ' trackdata = cluster_contour_tracker_update(videodata, template_struct, numIterations,...\n');
    fprintf(fid, '    coilIntensityCorrectionOn, coilSensitivityMatFileName, ...\n');
    fprintf(fid, '    newtonMethodOn, contourCleanUpOn, plotOn, notifyOn, frames, ...\n');
    fprintf(fid, '    parallelOn)\n');
    
    fprintf(fid, ' save(''../../%s_track.mat'',''trackdata'')\n',outFileName);
    %fprintf(fid, 'delete(gcp(''nocreate''))');
    fprintf(fid, 'catch\n');
    fprintf(fid, ' warning(''Experiment bypassed.'')\n'); 
    fprintf(fid, 'end\n');
    fclose(fid);
    
    %fprintf(fidout, 'cd %s/%s/experiment; main; clear\n',hpc_folder,pbsFolderCoreName);
    %fprintf(fidout2, 'cd %s/%s/experiment; matlab -nodesktop -nodisplay < main.m &> file.out &\n',hpc_folder,pbsFolderCoreName);
    %fprintf(fidout, 'cd %s/%s/experiment; main; clear\n',hpc_folder,pbsFolderCoreName);
    
    fprintf(fidout, 'cd %s/%s/experiment;\n',hpc_folder,pbsFolderCoreName);
    fprintf(fidout, 'c = parcluster(''24h_MultiNodeProfile'');\n');
    %fprintf(fidout, 'j = batch( c, ''main'', ''Pool'', 63, ''AttachedFiles'', {''gridkb.m'',''calckbkernel.m'',''kb.m'',''ift.m'',''interpft2.m''});\n');
    %fprintf(fidout, 'wait(j); delete(j);\n');
    fprintf(fidout, 'j = batch( c, ''main'',''Pool'', 63);\n');
    lineCount = lineCount + 1;
    
end;

%fprintf(fidout, 'quit\n');

fclose(fidin);
fclose(fidout);
%fclose(fidout2);



