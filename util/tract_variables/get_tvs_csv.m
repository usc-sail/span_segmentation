function get_tvs_csv(trackfile,csvfile)
%GET_TVS_CSV Summary of this function goes here
%   Detailed explanation goes here

contourdata = get_tvs_from_trackfile(trackfile);

fid = fopen(csvfile,'w');

n = size(contourdata.tv{1}.cd);

fprintf(fid,'Time,');
for i = 1:6
    
    fprintf(fid,'%s,',contourdata.tv{i}.name);
    
end;

fprintf(fid,'\n');

for i=1:n
    
    fprintf(fid, '%i,', round(1000*i/83.33));
    
    for j=1:6
        
        fprintf(fid,'%5.2f,',contourdata.tv{j}.cd(i));
        
    end;
    
    fprintf(fid,'\n');
    
end;
    
fclose(fid);

end

