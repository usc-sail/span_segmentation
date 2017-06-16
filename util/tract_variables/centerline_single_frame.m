function

load data.mat

frame=220;

fprintf('Trying to solve frame %i ...\n',frame);

vt_wall = get_vt_walls(contourdata, 1, frame);

for k1=0.55%0.70:-0.05:0.05
    
    
    for k2=0.70%0.95:-0.05:0.05
        
        
        success=1;
        
        try
            clf;
            iter = get_midline(vt_wall,k1,k2, 200, 5000,0);
        catch
            success = 0;
            fprintf('Failed frame %i using k1=%3.2f and k2=%3.2f\n',frame,k1,k2);
            continue;
        end;
        
        if success
            break;
        end;
        
        
    end;
    
    if success
        break;
    end;
end;

if success
    fprintf('SOLVED frame %i using k1=%3.2f and k2=%3.2f in %i iterations\n',frame,k1,k2, iter);
else
    fprintf('All attempts failed at frame %i\n',frame);
end;
