function [vt_wall] = get_vt_walls( contourdata, frameno )
%GET_VT_CONTOURS Get the internal and external vocal-tract walls of an MRI 
%   frame from a file in a contourdata structure
%
%   [internal_vt, external_vt] = GET_VT_CONTOURS(contourdata, fileno, frameno)
%
%   Asterios Toutios, Feb 25 2015

X=contourdata.X(frameno,:); 
Y=contourdata.Y(frameno,:); 
sectionsId=contourdata.SectionsID;

v=[X' Y'];

% 01 Epiglottis
% 02 Tongue
% 03 Incisor
% 04 Lower Lip
% 05 Jaw
% 06 Trachea
% 07 Pharynx
% 08 Upper Bound
% 09 Left Bound
% 10 Low Bound
% 11 Palate
% 12 Velum
% 13 Nasal Cavity
% 14 Nose
% 15 Upper Lip

ivt=[];
evt=[];

% Epiglottis
epiglottis=v(sectionsId==1,:);
[maxy_epiglottis,indmaxy_epiglottis]=max(epiglottis(:,2));
epiglottis=epiglottis(1:indmaxy_epiglottis,:);

% Tongue
tongue=v(sectionsId==2,:);
indtongue=find(tongue(:,2)>maxy_epiglottis);
tongue=tongue(indtongue,:);

% Lower Lip
lowerlip=v(sectionsId==4,:);
ivt=[epiglottis; tongue; lowerlip];

% Pharynx
pharynx=v(sectionsId==7,:);

% Hard Palate
palate=v(sectionsId==11,:);
palate=palate(end:-1:1,:);

% Velum
velum=v(sectionsId==12,:);
velum=velum(end:-1:1,:);
[miny_velum,indminy_velum]=min(velum(:,2));
velum=velum(indminy_velum:end,:);

% Upper Lip
upperlip=v(sectionsId==15,:);
upperlip=upperlip(end:-1:1,:);

lowpharynxpoint=InterX([-100,epiglottis(1,2);100,epiglottis(1,2)]',pharynx');
[~,indlpp]=min(lowpharynxpoint(1,:));
lpp=[lowpharynxpoint(1,indlpp),lowpharynxpoint(2,indlpp)];

for j=1:size(pharynx,1)-1
    
    if pharynx(j,1)<lpp(1) && pharynx(j+1,1)<lpp(1); continue; end;
    if pharynx(j,1)>lpp(1) && pharynx(j+1,1)>lpp(1); continue; end;
    if pharynx(j,2)<lpp(2) && pharynx(j+1,2)<lpp(2); continue; end;
    if pharynx(j,2)>lpp(2) && pharynx(j+1,2)>lpp(2); continue; end;
    
    pharynx_start=j+1; break;
    
end;

highpharynxpoint=InterX([-100,miny_velum;100,miny_velum]',pharynx');
[~,indhpp]=min(highpharynxpoint(1,:));
hpp=[highpharynxpoint(1,indlpp),highpharynxpoint(2,indlpp)];

for j=1:size(pharynx,1)-1
    
    if pharynx(j,1)<hpp(1) && pharynx(j+1,1)<hpp(1); continue; end;
    if pharynx(j,1)>hpp(1) && pharynx(j+1,1)>hpp(1); continue; end;
    if pharynx(j,2)<hpp(2) && pharynx(j+1,2)<hpp(2); continue; end;
    if pharynx(j,2)>hpp(2) && pharynx(j+1,2)>hpp(2); continue; end;
    
    pharynx_end=j; break;
    
end;

indpharynx=intersect(find(pharynx(:,2)>epiglottis(1,2)),find(pharynx(:,2)<miny_velum));
pharynx=[lpp; pharynx(pharynx_start:pharynx_end,:);hpp];

evt=[pharynx; velum; palate; upperlip];

%plot(ivt(:,1),ivt(:,2),'*-',evt(:,1),evt(:,2),'*-'); axis equal; axis off; shg;

vt_wall=struct('internal',ivt,'external',evt);

end

function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve
%   together with any self-intersection points.
%
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point.
%   Each factor of the 'C' arrays is essentially a matrix containing
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

%...Argument checks and assignment of L2
error(nargchk(1,2,nargin));
if nargin == 1,
    L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
else
    L2 = varargin{1}; hF = @le;
end

%...Preliminary stuff
x1  = L1(1,:)';  x2 = L2(1,:);
y1  = L1(2,:)';  y2 = L2(2,:);
dx1 = diff(x1); dy1 = diff(y1);
dx2 = diff(x2); dy2 = diff(y2);

%...Determine 'signed distances'
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

%...Obtain the segments where an intersection is expected
[i,j] = find(C1 & C2);
if isempty(i),P = zeros(2,0);return; end;

%...Transpose and prepare for output
i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

%...Solve system of eqs to get the common points
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
    dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';

    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end


