function struct_to_avi(videostruct,avifilename,startms,endms)
%STRUCT_TO_AVI Summary of this function goes here
%   Detailed explanation goes here

videodata = videostruct.frames;
rate = videostruct.framerate;

if nargin < 3
    startframe = 1;
    endframe = size(videodata,3);
else
    startframe = ceil(startms/1000*rate);
    endframe = floor(endms/1000*rate);
end;

writerObj = VideoWriter(avifilename);
writerObj.FrameRate = rate;
open(writerObj);

for movieframe=startframe:endframe
    
    img=videodata(:,:,movieframe);
    
    imagesc([-41.5 41.5], [-41.5 41.5],img); colormap(gray); axis equal; axis off; hold on;
    
    writeVideo(writerObj,getframe);
    
end;

close(writerObj);

end

