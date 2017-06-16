function [TRCD,x1,y1,x2,y2] = getPharyngealCD(Xphar,Yphar,Xroot,Yroot)
    [D,x1,y1,x2,y2] = spanPolyPolyDist_modified(Xroot,Yroot,Xphar,Yphar);
    TRCD = min(D);
    x1=x1(1); x2=x2(1); y1=y1(1); y2=y2(1);
end