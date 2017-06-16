function [X,Y,SectionsID,frames]=contour_to_table_resample(contourfile)

samplepoints{1}.i=[10 20 10 10 15 10];
samplepoints{2}.i=[30 20 10 5];
samplepoints{3}.i=[10 15 15 10 10];

load(contourfile)

ncontour = length(trackdata);

SectionsID=[];
X=[]; Y=[];
frames=[];

for i=1:ncontour
    
    if ~isempty(trackdata{i}.contours)

        segment=trackdata{i}.contours.segment;
        
        for s=1:(size(segment,2)-1)
            v = segment{s}.v;
            vnew=[];inew=[];
            sectionsId = segment{s}.i;
            for sId=1:max(sectionsId)
                ithis = find(sectionsId==sId);
                vthis = v(ithis,:);
                
                n=samplepoints{s}.i(sId);
                dv=diff(vthis,1,1);
                d=sqrt(sum(dv.^2,2));
                d=[0; cumsum(d)];
                di=linspace(0,d(end),n);
                vi=interp1(d,vthis,di,'pchip');
                
                vnew=[vnew; vi];
                inew=[inew; sId*ones(size(vi,1),1)];
            end;
            segment{s}.v = vnew;
            segment{s}.i = inew;
        end;
        
        for s=1
            section_id_row    = segment{s}.i;
            x_coord_row       = segment{s}.v(:,1);
            y_coord_row       = segment{s}.v(:,2);
            max_sec=max(section_id_row);
        end;
        
        for s=2
            section_id_row    = [section_id_row; max_sec+segment{s}.i];
            x_coord_row       = [x_coord_row; segment{s}.v(:,1)];
            y_coord_row       = [y_coord_row; segment{s}.v(:,2)];
            max_sec=max(section_id_row);
        end;
        
        for s=3
            section_id_row    = [section_id_row; max_sec+segment{s}.i];
            x_coord_row       = [x_coord_row; segment{s}.v(:,1)];
            y_coord_row       = [y_coord_row; segment{s}.v(:,2)];
        end;
        
        SectionsID=[SectionsID, section_id_row];
        X=[X,x_coord_row];
        Y=[Y,y_coord_row];
        frames=[frames, i];
              
    end;

end;
hold off;

SectionsID=SectionsID';
X=X';Y=Y';




    

