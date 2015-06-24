function points = find_old_grid(x_st,x_end,y_st,y_end,xEnlarged,yEnlarged)
    
    points = zeros(1,4);
    innerDomain = [x_st x_end y_st y_end];
    sx = length(xEnlarged);
    sy = length(yEnlarged);
    deviation = 10^(-11);
    
    for j=1:2
        
        gg = 0;
        for k = 1:sx
            if( innerDomain(j)-deviation<xEnlarged(k) && xEnlarged(k)<innerDomain(j) + deviation)
                gg = 1;
                points(j)=k;break;
            end
        end
        
        tt = 0;
        for k = 1:sy
            if( innerDomain(j+2) -deviation<yEnlarged(k) && yEnlarged(k)<innerDomain(2+j) + deviation)
                tt = 1;
                points(2+j)=k;break;
            end
        end
        
        if(gg == 0 || tt == 0)
                error('Cannot find point %d !',j);
        end
    end
end