function [dist_between_walls,dist_from_glottis] = get_sagittal_slice(centerpoints,touchpoints)

d=size(centerpoints,1);

dist_from_glottis = zeros(d,1);
dist_between_walls = zeros(d,1);

for j=1:d
    
    dist_between_walls(j) = sqrt((touchpoints(j).int(1) - touchpoints(j).ext(1)).^2 ...
        + (touchpoints(j).int(2) - touchpoints(j).ext(2)).^2);
    
end;

for j=2:d
    
    dist_from_glottis(j) = dist_from_glottis(j-1) ...
        + sqrt((centerpoints(j,1)-centerpoints(j-1,1))^2+(centerpoints(j,2)-centerpoints(j-1,2))^2);
    
end;