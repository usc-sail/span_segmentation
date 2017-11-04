function varargout = make_template(varargin)
% MAKE_TEMPLATE MATLAB GUI for drawing a vocal tract template for use
% with Erik Bresch's 2009 segmentation code
%
% Run as make_template(videodata, outputfile)
% The program will help you choose the frame in the .avi file that you
% will draw the template on.
% Two files {outputfile}.trk and {outputfile}.mat will be created
%
% Asterios Toutios 17-Sep-2014
% Last Modified by GUIDE v2.5 12-Dec-2016 10:57:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @make_template_OpeningFcn, ...
                   'gui_OutputFcn',  @make_template_OutputFcn, ...
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


% --- Executes just before make_template is made visible.
function make_template_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to make_template (see VARARGIN)

% Choose default command line output for make_template
handles.output = hObject;
handles.selected_curve=1;
handles.curves=struct('position',[]);
handles.position=[];
handles.h=[];
handles.template_struct_filename=varargin{1};


if ~isempty(varargin{3})
    coilSensitivityMatFileName=varargin{2};
    load(coilSensitivityMatFileName,'magnitudeCoilSensMap');
end;

% Update handles structure
guidata(hObject, handles);

videoData = varargin{1};
nFrames = size(videoData,3);

if ~isempty(varargin{3})
    videoData=videoData./repmat(magnitudeCoilSensMap,1,1,nFrames);
end;

%set(handles.slider1,'Max',nFrames);
%set(handles.slider1,'Min',1);

videoData=varargin{2}.frames;

handles.videoData=videoData;
handles.nFrames=size(handles.videoData,3);
handles.framerate=varargin{2}.framerate;

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

for i=1:15
    handles.curves(i).position=[];
end;

% Update handles structure
guidata(hObject, handles);

%load([handles.infilename,'.mat'],'template');
load(handles.template_struct_filename,'template_struct');
handles.itemplate=1;
template=template_struct(handles.itemplate).template;
handles.template_struct=template_struct;

% icurve=0;
% for isegment=1:3
%     for iv = 1:max(template.segment{1,isegment}.i);
%         icurve=icurve+1;
%         handles.curves(icurve).position=template.segment{1,isegment}.v(template.segment{1,isegment}.i==iv,:);
%         handles.curves(icurve).position=handles.imagewidth*handles.curves(icurve).position;
%         handles.curves(icurve).position(:,2)=-handles.curves(icurve).position(:,2)
%     end
% end;

for i=1:15
    handles.curves(i).position=template.curves(i).position;
end;

%image((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
colormap(gray); axis image; axis off;
hold on;

text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);

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


% UIWAIT makes make_template wait for user response (see UIRESUME)
% uiwait(handles.thisgui);


% --- Outputs from this function are returned to the command line.
function varargout = make_template_OutputFcn(hObject, eventdata, handles) 
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

text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);

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

text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);

for i=1:length(handles.curves)
    if ~isempty(handles.curves(i).position)
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
end;

hold off;

% --- Executes on button press in mri_save_template.
function mri_save_template_Callback(hObject, eventdata, handles)
% hObject    handle to mri_save_template (see GCBO)
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
text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);

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

frame=max(1,ceil(handles.nFrames*get(hObject,'Value')));
handles.frame=frame;
MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),MRI);
colormap(gray); axis image; axis off; hold on;
text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);
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



function frame_number_Callback(hObject, eventdata, handles)
% hObject    handle to frame_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_number as text
%        str2double(get(hObject,'String')) returns contents of frame_number as a double
frame=min(handles.nFrames,ceil(str2double(get(hObject,'String'))*handles.framerate/1000));
handles.frame=frame;
MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),MRI);
colormap(gray); axis image; axis off; hold on;
text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);
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
function frame_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in previous_frame.
function previous_frame_Callback(hObject, eventdata, handles)
% hObject    handle to previous_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame=max(handles.frame-1,1);
handles.frame=frame;
MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),MRI);
colormap(gray); axis image; axis off; hold on;
text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);
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


% --- Executes on button press in next_frame.
function next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame=min(handles.frame+1,handles.nFrames);
handles.frame=frame;
MRI=mat2gray(handles.videoData(:,:,frame),[handles.amin handles.amax]);
imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),MRI);
colormap(gray); axis image; axis off; hold on;
text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);
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



function template_number_Callback(hObject, eventdata, handles)
% hObject    handle to template_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of template_number as text
%        str2double(get(hObject,'String')) returns contents of template_number as a double
handles.itemplate=str2double(get(hObject,'String'));

template=handles.template_struct(handles.itemplate).template;

% icurve=0;
% for isegment=1:3
%     for iv = 1:max(template.segment{1,isegment}.i);
%         icurve=icurve+1;
%         handles.curves(icurve).position=template.segment{1,isegment}.v(template.segment{1,isegment}.i==iv,:);
%         handles.curves(icurve).position=handles.imagewidth*handles.curves(icurve).position;
%         handles.curves(icurve).position(:,2)=-handles.curves(icurve).position(:,2)
%     end
% end;

%image((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);

if ~isempty(template)
    imagesc(handles.MRI, 'XData', [-handles.axislimit handles.axislimit], 'YData', [-handles.axislimit handles.axislimit])
    colormap(gray); axis image; axis off;
    hold on;  
    for i=1:length(handles.curves)
        handles.curves(i)=template.curves(i);
        plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
    end;
    handles.template=template;
end;
text(-0.8*handles.axislimit,0.8*handles.axislimit,num2str(handles.frame),'Color','white','FontSize',24);

guidata(hObject, handles);

% imagesc((-handles.axislimit):(handles.axislimit), (-handles.axislimit):(handles.axislimit),handles.MRI);
% colormap(gray); axis image; axis off;
% hold on;
% 
% for i=1:length(handles.curves)
%     if ~isempty(handles.curves(i).position)
%         plot(handles.curves(i).position(:,1),handles.curves(i).position(:,2),'r','LineWidth',2);
%     end;
% end;

% --- Executes during object creation, after setting all properties.
function template_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to template_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
