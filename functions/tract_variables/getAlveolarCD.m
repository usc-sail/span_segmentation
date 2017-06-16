function [alveolarCD,x1,y1,x2,y2] = getAlveolarCD(Xpalate,Ypalate,Xtip,Ytip)
% Compute constriction degree at the alveolar place.
    [D,x1,y1,x2,y2] = spanPolyPolyDist_modified(Xtip,Ytip,Xpalate,Ypalate);
    alveolarCD = min(D);
    x1=x1(1); x2=x2(1); y1=y1(1); y2=y2(1);
end