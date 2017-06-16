function contourdata = get_contourdata_from_trackfile(trackfile)

contourdata=struct;
contourdata.X=[];
contourdata.Y=[];
contourdata.File={};
contourdata.SectionsID=[];
contourdata.Frames=[];

[X,Y,SectionsID,frames]=contour_to_table_resample(trackfile);

tmpX=[contourdata.X;X];
contourdata.X=tmpX;

tmpY=[contourdata.Y;Y];
contourdata.Y=tmpY;

N = size(X,1);
tmpFile = cell(1,N);

for i=1:N   
    tmpFile{i}=num2str(1);
end;

contourdata.File=[contourdata.File, tmpFile];

tmpFrames=[contourdata.Frames; frames'];
contourdata.Frames=tmpFrames;

contourdata.SectionsID=SectionsID(1,:);
