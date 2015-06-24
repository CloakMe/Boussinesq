function DrawEndSolution(x,y,h,zeroX,zeroY,al,bt,c,bigU,bigUTimeDer,compBox)
    if(x(1) == compBox.x_st)
        x_st = compBox.x_st;
        x_end = compBox.x_end;
        y_st = compBox.y_st;
        y_end = compBox.y_end;
    else
        x_st = compBox.x_st2;
        x_end = compBox.x_end2;
        y_st = compBox.y_st2;
        y_end = compBox.y_end2;   
    end
    
    mag = floor(10/h);
    xx=x(zeroX-mag+1:zeroX+mag-1); yy=y(zeroY-mag+1:zeroY+mag-1); 

    BigUNearCenter = bigU(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);
    bigUTimeDerNearCenter = bigUTimeDer(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);
    % bigICNearCenter = bigIC(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);

   
    
    figure(2);
    surf(xx,yy,(BigUNearCenter)');
    xlabel('x','FontSize', 14);
    ylabel('y', 'FontSize', 14);
    
    colormap(jet(512))
    %colorbar;
    %caxis([-.1, 2.4]);
    grid off;
    %axis equal;
    %axis off;
    view(31, 22);
    %view(0, 90);
    shading flat; 
    set(gca,'fontsize',14);
    
    
    figure(3);
    surf(xx,yy,( bigUTimeDerNearCenter'))%bbU'
    xlabel('x','FontSize', 14);
    ylabel('y', 'FontSize', 14);
    colormap(jet(512))
    %colorbar;
    grid off;
    view(31, 22);
    shading flat; 
    set(gca,'fontsize',14);
end
