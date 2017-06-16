function generate_pbs_batch(avifile, templateFileName, coilSensitivityFile, hpcFolder)

expSize=400;

[avipath,aviname]=fileparts(avifile);

pbsBaseFolderName = avipath;

mkdir(pbsBaseFolderName);

videostruct = avi_to_struct(avifile); videodataFull=videostruct.frames;

numberOfFrames=size(videodataFull,3);
numberOfBatches = ceil(numberOfFrames/expSize);

bashFileName = sprintf('%s/%s_cluster/script.sh', pbsBaseFolderName,aviname);

fidout = fopen(bashFileName,'w');

for ibatch = 1:numberOfBatches
   
    videoFileName = sprintf('%s_%02dof%02d',aviname,ibatch,numberOfBatches);
    
    frames=(expSize*(ibatch-1)+1):min(numberOfFrames,expSize*ibatch);
    
    videodata = videodataFull(:,:,frames);
    
    save([videoFileName,'.mat'],'videodata'); 
    
    pbsCoreFolderName = sprintf('%s_cluster/%s',aviname,videoFileName);
    
    generate_pbs(pbsBaseFolderName, pbsCoreFolderName, videoFileName,...
        1:length(frames), templateFileName, videoFileName,...
        coilSensitivityFile, hpcFolder);

    delete([videoFileName,'.mat']);
    
    fprintf(fidout, 'qsub -A lc_sn %s/pbs\n',videoFileName);  
    
end

fclose(fidout);



