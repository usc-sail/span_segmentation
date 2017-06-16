function [midpoint, cross_section, count] = get_midline(vt_wall, k1, k2, maxcircles, partiter, plotOn)
%GET_MIDLINE Get midline algorithm
%   Asterios Toutios, Feb 25 2015

% Load and upsample

ivt=InterpolateContourPoints2D(vt_wall.internal,100);
evt=InterpolateContourPoints2D(vt_wall.external,100);

%ivt=vt_wall.inetrnal;
%evt=vt_wall.external;

ivt=keep_vt_open(ivt,evt);

if plotOn
    plot(ivt(:,1),ivt(:,2),'-',evt(:,1),evt(:,2),'-'); axis equal; shg; hold on;%axis off; shg; hold on;
end;

circle = struct('center',[],'radius',[]);      % P1
circle_aux =struct('center',[],'radius',[]);    % P1'
circle_tmp = struct('center',[],'radius',[]);
cross_section = struct();
midpoint = zeros(100,2);

% First circle (P1)

circle(1).center = [evt(1,1) + ivt(1,1), evt(1,2) + ivt(1,2)]/2 ;

circle(1).radius = norm(evt(1,:) - circle(1).center,2); % euclidean distance

midpoint(1,:) = circle.center;

cross_section(1).int = ivt(1,:);
cross_section(1).ext = evt(1,:);

if plotOn
    plot(circle(1).center(1),circle(1).center(2),'r*');
    plotcircle(circle(1).center,circle(1).radius);
    pause(0.1);
end;

% First auxilliary circle (P1') is same as first P1

circle_aux(1).center = circle(1).center;
circle_aux(1).radius = k2*circle(1).radius;

% Iterative process

%k1=1;
%k2=0.65;
k3=0.01;
%k3=0.95;
k4=0.012;
%k4=1.05;
k5=0.09*pi/360;
k6=0.19*(pi/360);

% this will loop across the vocal tract


theta=0;

ii=2;

count = 0;

while(circle(ii-1).center(1) > (min(evt(:,1))+min(ivt(:,1)))/2)
%while(circle(ii-1).center(1) > -25)
    
    if ii > maxcircles
        error ('Max number of circles reached');
    end;
    
    circle_tmp.radius = k2*circle_aux(ii-1).radius;
    
    partcount = 0;

    while(1)
        
        circle_tmp.center = ...
            [circle_aux(ii-1).center(1) + circle_aux(ii-1).radius*cos(theta),...
            circle_aux(ii-1).center(2) + circle_aux(ii-1).radius*sin(theta)];
        
        [m1, tang_point1] = points_in_circle(circle_tmp, evt);
        [m2, tang_point2] = points_in_circle(circle_tmp, ivt);
        
        %[m1 m2]
        
        if (m1==1 && m2==1)
            break
        elseif (m1>1 && m2 >1)
            circle_tmp.radius = circle_tmp.radius - k3;
        elseif (m1==0 && m2==0);
            circle_tmp.radius = circle_tmp.radius + k4;
        elseif (m1 > m2)
            theta = theta+k5;
        else %(m1 < m2)
            theta = theta-k6;
        end
        
        count = count+1;
        partcount = partcount+1;
        
        if partcount > partiter
            error ('Max attempts for this circle reached');
        end;

        
    end;
    
    circle(ii) = circle_tmp;
    
    cross_section(ii).int = tang_point2;
    cross_section(ii).ext = tang_point1;
    
    midpoint(ii,:) = (tang_point1 + tang_point2)/2;
    
    if plotOn
        
        plotcircle(circle(ii).center,circle(ii).radius);
        plot(midpoint(ii,1),midpoint(ii,2),'k*');
        plot([cross_section(ii).int(1) cross_section(ii).ext(1)],[cross_section(ii).int(2) cross_section(ii).ext(2)],'k-');
        axis equal; shg; pause(0.1);
        
    end;
    
    % Create auxilliary (P') circle
    
    circle_aux(ii).radius = k2*circle(ii).radius;
    
    if ii>5
        anchor = circle(ii-5).center;
    else
        anchor=circle(1).center;
    end;
    
    
    [circle_aux(ii).center, theta] = get_aux_center ( circle(ii).center, midpoint(ii,:), k1*circle(ii).radius, anchor);
    
    if plotOn
        
        plotcircle_aux(circle_aux(ii).center,circle_aux(ii).radius);
        
    end;
    
    %axis equal; shg; pause(0.1);
    
    %circle_aux(ii)=circle(ii);
    
    ii=ii+1;
    
end;

midpoint = midpoint(1:(ii-1),:);

if plotOn
    
    plot(midpoint(:,1), midpoint(:,2),'k','LineWidth',2);
    axis equal; shg; pause(1);
    
    %     plotcircle(circle_tmp.center,circle_tmp.radius);
    %
    %     plot(midpoint(1,1),midpoint(1,2),'k*');
    %     plot([cross_section(1).int(1) cross_section(1).ext(1)],[cross_section(1).int(2) cross_section(1).ext(2)],'k-');
    
    
    hold off;
    
end;

%axis equal;

end

function plotcircle(center,radius)

t=linspace(0,2*pi,40);

plot(center(1)+radius*cos(t), center(2)+radius*sin(t), 'r');
plot(center(1), center(2), 'r*');

end

function plotcircle_aux(center,radius)

t=linspace(0,2*pi,40);

plot(center(1)+radius*cos(t), center(2)+radius*sin(t), 'g');
plot(center(1), center(2), 'g*');

end

function [m, tang_point] = points_in_circle(circle,line)

tmp = sqrt((line(:,1) - circle.center(1)).^2 + (line(:,2) - circle.center(2)).^2) < circle.radius ;

m = sum(tmp);
if m==1
    tang_point = line(tmp>0,:);
else
    tang_point = 0;
end;

end

function [aux_center, theta] = get_aux_center (prev_center, prev_midpoint, R, anchor)

x1 = prev_center(1);
y1 = prev_center(2);
x2 = prev_midpoint(1);
y2 = prev_midpoint(2);

a = (y1-y2)/(x1-x2);
b = (x1*y2-x2*y1)/(x1-x2);

x3a = x2 + R/sqrt(1+a^2);
y3a = a*x3a + b;

x3b = x2 - R/sqrt(1+a^2);
y3b = a*x3b + b;

da = norm([x3a y3a] - anchor,2);
db = norm([x3b y3b] - anchor,2);

if da > db
    aux_center = [x3a y3a];
    theta = atan(a);
else
    aux_center = [x3b y3b];
    theta = pi+atan(a);
end;


end

function ivtnew = keep_vt_open(ivt,evt)

ivtnew = ivt;

for i=round(size(ivt,1)/3):size(ivt,1)
    
    x_ivt = ivt(i,1);
    y_ivt = ivt(i,2);
    
    y_evt = interp1(evt(:,1),evt(:,2),x_ivt);
    
    if y_ivt > y_evt - 0.05
        
        y_ivt = y_evt - 0.05;
        
        ivtnew(i,2) = y_ivt;
        
    end
    
end

end
