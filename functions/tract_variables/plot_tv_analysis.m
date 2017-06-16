function plot_tv_analysis(trackfile,target_time)

contourdata = get_tvs_from_trackfile(trackfile);
fps = 83.33;
pixel = 0.24; %cm per pixel

n=size(contourdata.X,1);
time=1000*(1:n)/fps;

fig = figure('units','pixels','position',[100 100 1100 700])

if target_time>0;
    minframe = round(target_time*fps);
else
    minframe = find_narrowest_constriction(contourdata);
end;

[centerline,cross_sections] = centerline_from_frame(contourdata,minframe);
[slice_y,slice_x] = get_sagittal_slice(centerline,cross_sections);

subplot(2,3,1);
plot(time,contourdata.tv{1}.cd*pixel);
hold on; plot([minframe*1000/fps,minframe*1000/fps] ,[0 3],'r'); hold off;
title('bilabial constriction');
axis([0 time(end) 0 3]);
xlabel('Time (msec)');
ylabel('Constriction degree (cm)');

subplot(2,3,2);
plot(time,contourdata.tv{3}.cd*pixel);
hold on; plot([minframe*1000/fps,minframe*1000/fps] ,[0 3],'r'); hold off;
title('alveolar constriction');
axis([0 time(end) 0 3]);
xlabel('Time (msec)');
ylabel('Constriction degree (cm)');

subplot(2,3,4);
plot(time,contourdata.tv{4}.cd*pixel);
hold on; plot([minframe*1000/fps,minframe*1000/fps] ,[0 3],'r'); hold off;
title('palatal constriction');
axis([0 time(end) 0 3]);
xlabel('Time (msec)');
ylabel('Constriction degree (cm)');

subplot(2,3,5);
plot(time,contourdata.tv{5}.cd*pixel);
hold on; plot([minframe*1000/fps,minframe*1000/fps] ,[0 3],'r'); hold off;
title('velar constriction');
axis([0 time(end) 0 3]);
xlabel('Time (msec)');
ylabel('Constriction degree (cm)');

subplot(2,3,6)
plot(slice_x*pixel, slice_y*pixel,'m');
title('midsagittal slice at marked instance');
xlabel('Distance from glottis (cm)');
ylabel('Distance between walls (cm)');
axis([0 20 0 3]);

subplot(2,3,3);
hold on;
for section=1:15
    rrange = find(contourdata.SectionsID==section);
    plot(contourdata.X(minframe,rrange), contourdata.Y(minframe,rrange),'r-','LineWidth',2);
end;
title('vocal tract shape at constriction');
plot(centerline(:,1),centerline(:,2),'g');
for i=1:size(cross_sections,2)
    plot([cross_sections(i).int(1),cross_sections(i).ext(1)],...
        [cross_sections(i).int(2), cross_sections(i).ext(2)],'m');
end;
hold off;
axis([-40 20 -30 30]); axis off;

figfile = strrep(trackfile,'.mat','.png');
fig.PaperPositionMode = 'auto';

print(figfile,'-dpng','-r0');

end

function constriction_frame = find_narrowest_constriction(contourdata)

min_values=zeros(4,1);
min_indexes=zeros(4,1);

n=size(contourdata.tv{1}.cd);
rrange = round(n/4):round(3*n/4);

tv_indexes = [1 3 4 5]

for i=1:4
    
   [min_values(i),min_indexes(i)] = min(contourdata.tv{tv_indexes(i)}.cd(rrange));
    
end;

[~,constriction_location] = min(min_values);
constriction_frame = rrange(1)-1+min_indexes(constriction_location);

end

