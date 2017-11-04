function varargout = switch_make_template(varargin)
% SWITCH_MAKE_TEMPLATE MATLAB GUI for drawing a vocal tract template for use
% with Erik Bresch's 2009 segmentation code
%
% Run as switch_make_template(videodata, outputfile)
% The program will help you choose the frame in the .avi file that you
% will draw the template on.
% Two files {outputfile}.trk and {outputfile}.mat will be created
%
% Asterios Toutios 17-Sep-2014
% Last Modified by GUIDE v2.5 01-Aug-2016 13:02:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @switch_make_template_OpeningFcn, ...
                   'gui_OutputFcn',  @switch_make_template_OutputFcn, ...
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


% --- Executes just before switch_make_template is made visible.
function switch_make_template_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to switch_make_template (see VARARGIN)

% Choose default command line output for switch_make_template
handles.output = hObject;
handles.selected_curve=1;
handles.curves=struct('position',[]);
handles.position=[];
handles.h=[];
handles.infilename=varargin{2};
handles.outfilename=varargin{3};
handles.framerate=varargin{4};

if ~isempty(varargin{5})
    coilSensitivityMatFileName=varargin{5};
    load(coilSensitivityMatFileName,'magnitudeCoilSensMap');
end;

% Update handles structure
guidata(hObject, handles);

videoData = varargin{1};
nFrames = size(videoData,3);

if ~isempty(varargin{5})
    videoData=videoData./repmat(magnitudeCoilSensMap,1,1,nFrames);
end;

%set(handles.slider1,'Max',nFrames);
%set(handles.slider1,'Min',1);

handles.videoData=videoData;
handles.nFrames=nFrames;

amax=max(max(max(videoData)));
amin=min(min(min(videoData)));

handles.amax=amax;
handles.amin=amin;

handles.imagewidth = size(videoData,1);
handles.axislimit = ((handles.imagewidth-1)/2);

MRI=mat2gray(videoData(:,:,1),[amin amax]);
handles.MRI=MRI;

%image((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;

for i=1:15
    handles.curves(i).position=[];
end;

% Update handles structure
guidata(hObject, handles);

load([handles.infilename,'.mat'],'template');

for i=1:15
    handles.curves(i).position=template.curves(i).position;
end;

%image((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
hold on;

for i=1:length(handles.curves)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
end;

guidata(hObject, handles);

%h=impoly;
%position=wait(h);
%
%handles.h=h;


% UIWAIT makes switch_make_template wait for user response (see UIRESUME)
% uiwait(handles.thisgui);


% --- Outputs from this function are returned to the command line.
function varargout = switch_make_template_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on button press in mri_draw.
function mri_draw_Callback(hObject, eventdata, handles)
% hObject    handle to mri_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
delete(handles.h);
handles.h=impoly('Closed',false);
%handles.position=getPosition(handles.h);
%handles.position=position;
guidata(hObject, handles);

% --- Executes on button press in mri_edit.
function mri_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mri_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected_curve=handles.selected_curve;
if ~isempty(handles.curves(selected_curve).position)
    initposition=handles.curves(selected_curve).position;
else
    initposition=[];
end;
%handles.curves(selected_curve).position=position;

imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
colormap(gray); axis image; axis off;
hold on;

for i=setdiff(1:length(handles.curves),selected_curve)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
end;


handles.h=impoly(gca,initposition,'Closed',false);
%position2=getPosition(h);
%handles.position=position2;
guidata(hObject, handles);

%hold off;



% --- Executes on button press in mri_save.
function mri_save_Callback(hObject, eventdata, handles)
% hObject    handle to mri_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.h)
    handles.position=getPosition(handles.h);
    selected_curve=handles.selected_curve;
    %position=handles.position;
    handles.curves(selected_curve).position=handles.position;
    guidata(hObject, handles);
end;

imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
colormap(gray); axis image; axis off;
hold on;

for i=1:length(handles.curves)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
end;

hold off;

% --- Executes on button press in mri_save_all_exit.
function mri_save_all_exit_Callback(hObject, eventdata, handles)
% hObject    handle to mri_save_all_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fid = fopen([handles.outfilename,'.trk'],'w');
for i=1:15
    
    if ~isempty(handles.curves(i).position)
        
        ctrx=handles.curves(i).position(:,1);
        ctry=handles.curves(i).position(:,2);
        
        fprintf(fid,'im1\n');
        fprintf(fid,'c%i\n\n',i);
        for j=1:length(ctrx)
            fprintf(fid,'%6.3f,%6.3f\n',ctrx(j),-ctry(j));
        end;
        fprintf(fid,'\n');
        
    end;
    
end;

fclose(fid);

template=struct('curves',[],'MRI',[]);

template.curves = handles.curves;
template.MRI = handles.MRI;

save([handles.outfilename,'.mat'],'template');

switch_convert_template([handles.outfilename,'_converted.mat'],[handles.outfilename,'.trk'],1,handles.imagewidth,133,1);

guidata(hObject, handles);

close(handles.thisgui);


% --- Executes on selection change in curve_selection.
function curve_selection_Callback(hObject, eventdata, handles)
% hObject    handle to curve_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns curve_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from curve_selection
handles.selected_curve=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function curve_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curve_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected_curve=handles.selected_curve;
handles.curves(selected_curve).position=[];
guidata(hObject, handles);

imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
colormap(gray); axis image; axis off;
hold on;

for i=1:length(handles.curves)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
end;

hold off;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

frame=ceil(handles.nFrames*get(hObject,'Value'));
MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),MRI);
colormap(gray); axis image; axis off; hold on;
handles.MRI=MRI;
% for i=1:15
%     handles.curves(i).position=[];
% end;
for i=1:length(handles.curves)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
frame=min(handles.nFrames,ceil(str2double(get(hObject,'String'))*handles.framerate/1000))
MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),MRI);
colormap(gray); axis image; axis off; hold on;
handles.MRI=MRI;
% for i=1:15
%     handles.curves(i).position=[];
% end;
for i=1:length(handles.curves)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
end;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
