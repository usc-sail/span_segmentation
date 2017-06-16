function contourdata = get_tvs_from_trackfile(trackfile)

ds = 1/6; % Make tissue-air boundary outlines denser by a factor of six.

contourdata = get_contourdata_from_trackfile(trackfile);
n=size(contourdata.X,1);

LA=zeros(n,1); ulx=LA; uly=LA; llx=LA; lly=LA;
VEL=LA; velumx1=LA; velumy1=LA; pharynxx1=LA; pharynxy1=LA;
alveolarCD=LA; alveolarx=LA; alveolary=LA; tonguex2=LA; tonguey2=LA;
palatalCD=LA; palatalx=LA; palataly=LA; tonguex3=LA; tonguey3=LA;
velarCD=LA; velumx2=LA; velumy2=LA; tonguex1=LA; tonguey1=LA;
pharyngealCD=LA; pharynxx2=LA; pharynxy2=LA; tonguex4=LA; tonguey4=LA; 

files = unique(contourdata.File);
nFiles = length(files);
k=1;
for i=1:nFiles
    for j=contourdata.Frames(strcmp(contourdata.File,files(i)))'
        
        [X,Y,Xul,Yul,Xll,Yll,Xtongue,Ytongue,Xalveolar,...
            Yalveolar,Xpalatal,Ypalatal,Xvelum,Yvelum,...
            Xvelar,Yvelar,Xphar,Yphar,Xepig,Yepig] = vtSeg(contourdata,...
            files(i),j,ds);
        
        [LA(k),ulx(k),uly(k),llx(k),lly(k)] = getLA(Xul,Yul,Xll,Yll);
        [VEL(k),velumx1(k),velumy1(k),pharynxx1(k),pharynxy1(k)] = getVEL(Xvelum,Yvelum,Xphar,Yphar);
        [alveolarCD(k),alveolarx(k),alveolary(k),tonguex2(k),tonguey2(k)] = getAlveolarCD(Xalveolar,Yalveolar,Xtongue,Ytongue);
        [palatalCD(k),palatalx(k),palataly(k),tonguex3(k),tonguey3(k)] = getPalatalCD(Xpalatal,Ypalatal,[Xtongue, Xepig],[Ytongue,Yepig]);
        [velarCD(k),velumx2(k),velumy2(k),tonguex1(k),tonguey1(k)] = getVelarCD(Xvelar,Yvelar,[Xtongue,Xepig],[Ytongue,Yepig]);
        [pharyngealCD(k),pharynxx2(k),pharynxy2(k),tonguex4(k),tonguey4(k)] = getPharyngealCD(Xphar,Yphar,[Xtongue,Xepig],[Ytongue,Yepig]);
        
        k=k+1;
    end
end

contourdata.tv{1}.cd=LA; contourdata.tv{1}.in=[llx lly]; contourdata.tv{1}.out=[ulx uly];
contourdata.tv{2}.cd=VEL; contourdata.tv{2}.in=[velumx1 velumy1]; contourdata.tv{2}.out=[pharynxx1 pharynxy1];
contourdata.tv{3}.cd=alveolarCD; contourdata.tv{3}.in=[tonguex2 tonguey2]; contourdata.tv{3}.out=[alveolarx alveolary];
contourdata.tv{4}.cd=palatalCD; contourdata.tv{4}.in=[tonguex3 tonguey3]; contourdata.tv{4}.out=[palatalx palataly];
contourdata.tv{5}.cd=velarCD; contourdata.tv{5}.in=[tonguex1 tonguey1]; contourdata.tv{5}.out=[velumx2 velumy2];
contourdata.tv{6}.cd=pharyngealCD; contourdata.tv{6}.in=[tonguex4 tonguey4]; contourdata.tv{6}.out=[pharynxx2 pharynxy2];

contourdata.tv{1}.name = 'Bilabial';
contourdata.tv{2}.name = 'Velophrayngeal';
contourdata.tv{3}.name = 'Alveolar';
contourdata.tv{4}.name = 'Palatal';
contourdata.tv{5}.name = 'Velar';
contourdata.tv{6}.name = 'Pharyngeal';