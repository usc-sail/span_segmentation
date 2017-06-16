function get_sagittal_csv(trackfile,csvfile,target_time)
%GET_TVS_CSV Summary of this function goes here
%   Detailed explanation goes here

fps = 83.33;

contourdata = get_tvs_from_trackfile(trackfile);

if target_time>0;
    minframe = round(target_time*fps);
else
    minframe = find_narrowest_constriction(contourdata);
end;

[centerline,cross_sections] = centerline_from_frame(contourdata,minframe);
pixel = 0.24;
[slice_y,slice_x] = get_sagittal_slice(centerline,cross_sections);

plot(slice_x,slice_y)

fid = fopen(csvfile,'w');

n = length(slice_x);

for i=1:n
        
        fprintf(fid,'%5.2f, %5.2f\n',pixel*slice_x(i),pixel*slice_y(i));
    
end;
    
fclose(fid);

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