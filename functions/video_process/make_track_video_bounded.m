function make_track_video_bounded(avifile,contourfile)

videostruct = avi_to_struct(avifile);
videodata = videostruct.frames;
rate = videostruct.framerate;

load(contourfile,'trackdata');

close all;

set(gca,'NextPlot','replaceChildren');

nmovie = size(videodata,3);
width = size(videodata,1);

trackvideofile = strrep(contourfile,'.mat','.avi');
writerObj = VideoWriter(trackvideofile);
writerObj.FrameRate = videostruct.framerate;
open(writerObj);

%readObj = VideoReader(videofile);

dataindex=1;

close all; figure;

for movieframe=1:nmovie
    
    if dataindex <= size(trackdata,2)
        dataframe=trackdata{dataindex}.frameNo;
    else
        dataframe=0;
    end;
    

    if ~isempty(dataframe)
        
        img=videodata(:,:,movieframe);
     
        if movieframe==dataframe
        imshow(mat2gray(img), 'XData',[-(width-1)/2 (width-1)/2],'YData',[-(width-1)/2 (width-1)/2],'Border','tight','InitialMagnification',500); hold on;
            
            segment=trackdata{dataindex}.contours.segment;
            
            for s=1:(size(segment,2)-1)
                sectionsId = segment{s}.i;
                v          = segment{s}.v;
                colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
                for sId=1:max(sectionsId)
                    plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),colors(sId),'LineWidth',4); hold on;
                end;
            end;
            dataindex=dataindex+1;
               
        hold off;
        
        drawnow;
        writeVideo(writerObj,getframe);
            
        end;

        
    end;
    
    
end;

close(writerObj);

wavfile = strrep(avifile,'.avi','.wav');
[y,fs] = audioread(wavfile);

firstframe = trackdata{1}.frameNo;
lastframe = trackdata{dataindex-1}.frameNo;

first_ms = firstframe/rate;
last_ms = lastframe/rate;

first_sample = round(fs*first_ms);
last_sample = round(fs*last_ms);

y_cut = y(first_sample:last_sample);

outwavfile = strrep(contourfile,'.mat','.wav');

audiowrite(outwavfile,y_cut,fs);

outVideo = strrep(contourfile,'.mat','_with_audio.avi');

s = sprintf('ffmpeg -i %s -i %s -c:v copy -c:a aac -pix_fmt yuv420p -strict experimental %s', ...
    trackvideofile, outwavfile, outVideo);

system(s);
