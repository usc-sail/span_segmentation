function varargout = inspect_tracks2(varargin)
% INSPECT_TRACKS2 MATLAB GUI for drawing a vocal tract template for use
% with Erik Bresch's 2009 segmentation code
%
% Run as inspect_tracks2(videodata, outputfile)
% The program will help you choose the frame in the .avi file that you
% will draw the template on.
% Two files {outputfile}.trk and {outputfile}.mat will be created
%
% Asterios Toutios 17-Sep-2014
% Last Modified by GUIDE v2.5 16-Jan-2017 16:06:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inspect_tracks2_OpeningFcn, ...
                   'gui_OutputFcn',  @inspect_tracks2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before inspect_tracks2 is made visible.
function inspect_tracks2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inspect_tracks2 (see VARARGIN)

% Choose default command line output for inspect_tracks2
handles.output = hObject;
handles.selected_curve=1;
handles.curves=struct('position',[]);
handles.position=[];
handles.h=[];
handles.template_struct_filename=varargin{1};
handles.hide = 0;


if ~isempty(varargin{3})
    coilSensitivityMatFileName=varargin{2};
    load(coilSensitivityMatFileName,'magnitudeCoilSensMap');
end;

% Update handles structure
guidata(hObject, handles);


videoData=varargin{1}.frames;
%videoData = varargin{1};
trackdata = varargin{2};
nFrames = size(videoData,3);

if ~isempty(varargin{3})
    videoData=videoData./repmat(magnitudeCoilSensMap,1,1,nFrames);
end;

%set(handles.slider1,'Max',nFrames);
%set(handles.slider1,'Min',1);


handles.videoData=videoData;
handles.trackdata=trackdata;
handles.nFrames=size(handles.videoData,3);
handles.framerate=varargin{1}.framerate;

nTracks = size(trackdata,2);
handles.trackframes = zeros(nTracks,1);

for i=1:nTracks
    handles.trackframes(i)=trackdata{i}.frameNo;
end;
    
    
amax=max(max(max(videoData)));
amin=min(min(min(videoData)));

handles.amax=amax;
handles.amin=amin;

handles.imagewidth = size(videoData,1);
handles.axislimit = ((handles.imagewidth-1)/2);

MRI=mat2gray(videoData(:,:,1),[amin amax]);
handles.MRI=MRI;
handles.frame=1;

imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
axis
handles.axislimit

%segment=handles.trackdata{handles.frame}.contours.segment;



% Update handles structure
%guidata(hObject, handles);

%load([handles.infilename,'.mat'],'template');
%load(handles.template_struct_filename,'template_struct');
%handles.itemplate=1;
%template=template_struct(handles.itemplate).template;
%handles.template_struct=template_struct;
%handles.numberOfTemplates=size(template_struct,2);

%refresh_plot(handles)

imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
hold on;


trackframe = find(handles.trackframes == handles.frame);
s=sprintf('Frame %i \n%i ms',handles.frame, round(handles.frame/handles.framerate*1000));
    text(-0.9*handles.axislimit,0.8*handles.axislimit,s,'Color','white','FontSize',24);

if ~isempty(trackframe) && (~handles.hide)
    segment=handles.trackdata{trackframe}.contours.segment;
    
    s = sprintf ('Template %i', handles.trackdata{trackframe}.template);
    text(-0.5*handles.axislimit,0.8*handles.axislimit,s,'Color','red','FontSize',24);
    
    for s=1:(size(segment,2)-1)
        sectionsId = segment{s}.i;
        v          = segment{s}.v;
        colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
        for sId=1:max(sectionsId)
            plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),colors(sId),'LineWidth',4); hold on;
        end;
    end;
    
end;

guidata(hObject, handles);


% UIWAIT makes inspect_tracks2 wait for user response (see UIRESUME)
% uiwait(handles.thisgui);


% --- Outputs from this function are returned to the command line.
function varargout = inspect_tracks2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on button press in previous_frame.
function previous_frame_Callback(hObject, eventdata, handles)
% hObject    handle to previous_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame=max(handles.frame-1,1);
handles.frame=frame;

handles.MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
hold on;


trackframe = find(handles.trackframes == handles.frame);
s=sprintf('Frame %i \n%i ms',handles.frame, round(handles.frame/handles.framerate*1000));
    text(-0.9*handles.axislimit,0.8*handles.axislimit,s,'Color','white','FontSize',24);

if ~isempty(trackframe) && (~handles.hide)
    segment=handles.trackdata{trackframe}.contours.segment;
    
    s = sprintf ('Template %i', handles.trackdata{trackframe}.template);
    text(-0.5*handles.axislimit,0.8*handles.axislimit,s,'Color','red','FontSize',24);
    
    for s=1:(size(segment,2)-1)
        sectionsId = segment{s}.i;
        v          = segment{s}.v;
        colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
        for sId=1:max(sectionsId)
            plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),colors(sId),'LineWidth',4); hold on;
        end;
    end;
    
end;

guidata(hObject, handles);


% --- Executes on button press in next_frame.
function next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame=min(handles.frame+1,handles.nFrames);
handles.frame=frame;

handles.MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
hold on;


trackframe = find(handles.trackframes == handles.frame);
s=sprintf('Frame %i \n%i ms',handles.frame, round(handles.frame/handles.framerate*1000));
    text(-0.9*handles.axislimit,0.8*handles.axislimit,s,'Color','white','FontSize',24);

if ~isempty(trackframe) && (~handles.hide)
    segment=handles.trackdata{trackframe}.contours.segment;
    
    s = sprintf ('Template %i', handles.trackdata{trackframe}.template);
    text(-0.5*handles.axislimit,0.8*handles.axislimit,s,'Color','red','FontSize',24);
    
    for s=1:(size(segment,2)-1)
        sectionsId = segment{s}.i;
        v          = segment{s}.v;
        colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
        for sId=1:max(sectionsId)
            plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),colors(sId),'LineWidth',4); hold on;
        end;
    end;
    
end;

guidata(hObject, handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

frame=max(1,ceil(handles.nFrames*get(hObject,'Value')));
handles.frame=frame;

handles.MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
hold on;


trackframe = find(handles.trackframes == handles.frame);
s=sprintf('Frame %i \n%i ms',handles.frame, round(handles.frame/handles.framerate*1000));
    text(-0.9*handles.axislimit,0.8*handles.axislimit,s,'Color','white','FontSize',24);

if ~isempty(trackframe) && (~handles.hide)
    segment=handles.trackdata{trackframe}.contours.segment;
    
    s = sprintf ('Template %i', handles.trackdata{trackframe}.template);
    text(-0.5*handles.axislimit,0.8*handles.axislimit,s,'Color','red','FontSize',24);
    
    for s=1:(size(segment,2)-1)
        sectionsId = segment{s}.i;
        v          = segment{s}.v;
        colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
        for sId=1:max(sectionsId)
            plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),colors(sId),'LineWidth',4); hold on;
        end;
    end;
    
end;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in hidetracks.
function hidetracks_Callback(hObject, eventdata, handles)
% hObject    handle to hidetracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hidetracks
handles.hide = get(hObject,'Value');

frame=handles.frame;

handles.MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;

trackframe = find(handles.trackframes == handles.frame);
s=sprintf('Frame %i \n%i ms, template %i',handles.frame, round(handles.frame/handles.framerate*1000));
    text(-0.9*handles.axislimit,0.8*handles.axislimit,s,'Color','white','FontSize',24);

if ~isempty(trackframe) && (~handles.hide)
    segment=handles.trackdata{trackframe}.contours.segment;
    
    s = sprintf ('Template %i', handles.trackdata{trackframe}.template);
    text(-0.5*handles.axislimit,0.8*handles.axislimit,s,'Color','red','FontSize',24);
    
    for s=1:(size(segment,2)-1)
        sectionsId = segment{s}.i;
        v          = segment{s}.v;
        colors = ['r' 'g' 'b' 'y' 'c' 'm' 'k'];
        for sId=1:max(sectionsId)
            plot( v(sectionsId==sId,1),-v(sectionsId==sId,2),colors(sId),'LineWidth',4); hold on;
        end;
    end;
    
end;

guidata(hObject, handles);