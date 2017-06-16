function [palatalCD,x1,y1,x2,y2] = getPalatalCD(Xpalatal,Ypalatal,Xbody,Ybody)
% Compute constriction degree at the palatal place.
    [D,x1,y1,x2,y2] = spanPolyPolyDist_modified(Xbody,Ybody,Xpalatal,Ypalatal);
    palatalCD = min(D); 
    x1=x1(1); x2=x2(1); y1=y1(1); y2=y2(1);
end