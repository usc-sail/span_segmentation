function generate_pbs_from_file(pbsBaseFolderName, scriptName, segmentFile, firstLine, lastLine, frameRate, templateFileName, coilSensitivityFile, hpcFolder, account)
lineCount = 0;

bashFileName = sprintf('%s/%s.sh', pbsBaseFolderName, scriptName);
 
fid = fopen(segmentFile);
fidout = fopen(bashFileName,'w');
fprintf(fidout, '#!/bin/bash\n');

while lineCount<firstLine - 1
    
    fgetl(fid);
    lineCount = lineCount + 1;
    
end;

while lineCount<lastLine
    
    str = fgetl(fid);
    celldata = textscan(str,'%s','Delimiter',',');
    strdata = celldata{1};
    videoFileName = char(strdata(1));
    pbsCoreFolderName = char(strdata(2));
    outFileName = pbsCoreFolderName;
      
    frames=[];
    
    for segment = 1:(length(strdata)-2)/2
    
    segmentstartframe = floor(frameRate/1000*str2double(strdata(2*segment+1)));
    segmentendframe = ceil(frameRate/1000*str2double(strdata(2*segment + 2)));
    
    frames = [frames, segmentstartframe:segmentendframe];
    
    %avifile
    %videostruct = avi_to_struct(avifile); videodata=videostruct.frames;
    %videoFileName='tmp';
    %save([videoFileName,'.mat'],'videodata');
    
    end;
    
    generate_pbs(pbsBaseFolderName, pbsCoreFolderName, videoFileName,...
        frames, templateFileName, outFileName,...
        coilSensitivityFile, hpcFolder);
    
    fprintf(fidout, 'qsub -A %s %s/pbs\n',account,pbsCoreFolderName);
    
    lineCount = lineCount + 1;
    
end;

fclose(fid);
fclose(fidout);



