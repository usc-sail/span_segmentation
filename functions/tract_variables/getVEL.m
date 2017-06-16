function [VEL,x1,y1,x2,y2] = getVEL(Xvelum,Yvelum,Xphar,Yphar)
% Compute velic aperture.
    [D,x1,y1,x2,y2] = spanPolyPolyDist_modified(Xvelum,Yvelum,Xphar,Yphar);
    VEL = min(D);
    x1=x1(1); x2=x2(1); y1=y1(1); y2=y2(1);
end