function videostruct = avi_to_struct(avifile)

% function videostruct = avi_to_struct(avifile)
%
% Converts an avi file to a video structure
%
% Asterios Toutios, Sep 9 2015

% Determine number of frames

nmovie=0;

vidObj = VideoReader(avifile);

try
    
    nmovie = get(vidObj, 'NumberOfFrames');
    
catch
    
    while hasFrame(vidObj)
        readFrame(vidObj);
        nmovie = nmovie + 1;
    end
    
end

% Determine dimensions of frame an Framerate

xdim = vidObj.Height;
ydim = vidObj.Width;

framerate = vidObj.FrameRate;

% Reset

vidObj = VideoReader(avifile);

% Initialize

videostruct.frames = zeros(xdim,ydim,nmovie);
videostruct.framerate = framerate;

% Read

try
    
    for i=1:nmovie
        
        f=readFrame(vidObj);
        
        if ismatrix(f)
            videostruct.frames(:,:,i) = f;
        else
            videostruct.frames(:,:,i) = rgb2gray(f);
        end;
        
    end;
    
catch
    
    f=read(vidObj);
    for i=1:nmovie
        
        frame = f(:,:,:,i);
        
        if ismatrix(frame)
            videostruct.frames(:,:,i) = frame;
        else
            videostruct.frames(:,:,i) = rgb2gray(frame);
        end;
        
    end;

end;
    


end

