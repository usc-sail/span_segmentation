function cluster_generate_scripts_split(localFolderName, scriptName, segmentFile, firstLine, lastLine, templateFileName, coilSensitivityFile, hpc_folder, funFolder)
lineCount = 0;
frames_per_split = 250;
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
    
    videostruct=avi_to_struct(videoFileName); videodata=videostruct.frames;
    nframes = size(videodata,3);
    nsplits = ceil(nframes/frames_per_split);
    
    for split = 1:nsplits
        
        minframe = (split-1)*frames_per_split + 1;
        maxframe = min(split*frames_per_split,nframes);
        
        
        pbsFolderName = sprintf('%s/%s_%i_of_%i',localFolderName,pbsFolderCoreName,split,nsplits);
        
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
        
        
        save([pbsFolderName,'/experiment/videodata.mat'],'videodata');
        
        fid=fopen([pbsFolderName,'/experiment/main.m'],'w');
        fprintf(fid, 'warning(''off'',''all'');\n');
        fprintf(fid, 'load template.mat;\n');
        fprintf(fid, 'load videodata.mat;\n');
        fprintf(fid, 'frames = %i:%i;\n',minframe,maxframe);
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
        fprintf(fid, 'coilSensitivityMatFileName = ''coilSensitivity.mat'';\n');
        
        fprintf(fid, 'try\n');
        fprintf(fid, ' trackdata = cluster_contour_tracker_update(videodata, template_struct, numIterations,...\n');
        fprintf(fid, '    coilIntensityCorrectionOn, coilSensitivityMatFileName, ...\n');
        fprintf(fid, '    newtonMethodOn, contourCleanUpOn, plotOn, notifyOn, frames, ...\n');
        fprintf(fid, '    parallelOn)\n');
        
        fprintf(fid, ' save(''../../%s_%i_of_%i_track.mat'',''trackdata'')\n',pbsFolderCoreName,split,nsplits);
        fprintf(fid, 'catch\n');
        fprintf(fid, ' warning(''Experiment bypassed.'')\n');
        fprintf(fid, 'end\n');
        fclose(fid);
        
        fprintf(fidout, 'cd %s/%s_%i_of_%i/experiment;\n',hpc_folder,pbsFolderCoreName,split,nsplits);
        fprintf(fidout, 'c = parcluster(''24h_MultiNodeProfile'');\n');
        fprintf(fidout, 'j = batch( c, ''main'',''Pool'', 63);\n');
        
        
    end;
    
    lineCount = lineCount + 1;
    
end;

fclose(fidin);
fclose(fidout);



