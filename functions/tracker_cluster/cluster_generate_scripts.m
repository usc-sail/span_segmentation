function cluster_generate_scripts(localFolderName, scriptName, segmentFile, firstLine, lastLine, frameRate, templateFileName, coilSensitivityFile, hpc_folder, funFolder)
lineCount = 0;
frames =[];
scriptFileName = sprintf('%s/%s.m', localFolderName, scriptName);

fidin = fopen(segmentFile);
fidout = fopen(scriptFileName,'w');

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
    copyfile([funFolder,'/cluster_contour_tracker.m'],[pbsFolderName,'/experiment']);
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
    fprintf(fid, 'workers = 64;\n');
    fprintf(fid, 'coilSensitivityMatFileName = ''coilSensitivity.mat'';\n');
    
    fprintf(fid, 'trackdata = cluster_contour_tracker(videodata, template_struct, numIterations,...\n');
    fprintf(fid, '    coilIntensityCorrectionOn, coilSensitivityMatFileName, ...\n');
    fprintf(fid, '    newtonMethodOn, contourCleanUpOn, plotOn, notifyOn, frames, ...\n');
    fprintf(fid, '    parallelOn, workers)\n');
    
    fprintf(fid, 'save(''../../%s_track.mat'',''trackdata'')\n',outFileName);
    %fprintf(fid, 'delete(gcp(''nocreate''))');
    
    fclose(fid);
    
    fprintf(fidout, 'cd %s/%s/experiment; main; clear\n',hpc_folder,pbsFolderCoreName);
    
    lineCount = lineCount + 1;
    
end;

fprintf(fidout, 'quit\n');

fclose(fidin);
fclose(fidout);



