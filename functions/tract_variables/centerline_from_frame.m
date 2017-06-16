function [centerline,cross_sections] = centerline_from_frame(contourdata,frame)

for ii = frame
    
    vt_wall = get_vt_walls(contourdata, ii);
    
    for k1=0.55:0.05:1.25
       
        for k2=0.6:0.05:0.95
          
            success=1;
            
            try
                clf;
                [midpoint, cross_section, iter] = get_midline(vt_wall,k1,k2, 200, 5000,1);
            catch
                success = 0;
                fprintf('Failed frame using k1=%3.2f and k2=%3.2f\n',k1,k2);
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
        fprintf('SOLVED frame using k1=%3.2f and k2=%3.2f in %i iterations\n',k1,k2, iter);
        centerline = midpoint;
        cross_sections = cross_section;
    else
        fprintf('All attempts failed at ths frame');
    end;
    
end;