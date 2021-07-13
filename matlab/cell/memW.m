function [content_w,MolarDen_w]=memW(a)
    Den_mem=2000;
    M_mem=1.1;
    if 0<a<=1
        content_w=0.043+17.81*a-39.85*a^2+36*a^3;
    else if a<=3
        content_w=14+1.4*(a-1);
        else
        content_w=16.8;
        end
    end
    MolarDen_w=Den_mem*content_w/M_mem;
 
end