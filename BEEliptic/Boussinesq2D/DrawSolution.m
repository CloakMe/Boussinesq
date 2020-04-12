function DrawSolution(x,y,h,zeroX,zeroY,al,bt,c,theta,bigU,bigUTimeDer,bigIC,U,compBox,secondDerivative)
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
    
    mag = max(20, floor(5/h));
    xx=x(zeroX-mag+1:zeroX+mag-1); yy=y(zeroY-mag+1:zeroY+mag-1); 

    BigUNearCenter = bigU(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);
    bigUTimeDerNearCenter = bigUTimeDer(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);
    % bigICNearCenter = bigIC(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);


    figure(2);
    mesh(xx,yy,( bigUTimeDerNearCenter'))%bbU'
    xlabel('x');    ylabel('y');
    title('time derivative')

    figure(3);
    mesh(xx,yy,(BigUNearCenter)');
    xlabel('x');    ylabel('y');
    title('end solution')

    Q = 41;
    viewTypeX = 90;
    viewTypeY = 90;
    figure(21)
    mesh(x, y(1:Q), bigU(:,1:Q)');
    view( viewTypeX, viewTypeY );
    title('Bottom domain boundary near y_{start}');
    xlabel('x');            ylabel('y');
    figure(22)
    mesh(x(1:Q), y, bigU(1:Q,:)');
    view( viewTypeX, viewTypeY );    
    title('Left domain boundary near x_{start}');
    xlabel('x');            ylabel('y');

    figure(23)
    mesh(x, y(end-Q:end), bigU(:,end-Q:end)');
    view( viewTypeX, viewTypeY );    
    title('Top domain boundary near y_{end}');
    xlabel('x');            ylabel('y');    
    figure(24)
    mesh(x(end-Q:end), y, bigU(end-Q:end,:)');
    view( viewTypeX, viewTypeY );
    title('Right domain boundary near x_{end}');
    xlabel('x');            ylabel('y');
    
    Yk = bigU;
    bigZeroMatrix=zeros(size(Yk));
    dyy_yk = bt*c^2*YDer(Yk,secondDerivative);
    deltaU = Delta(Yk,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,secondDerivative);
    residual =  Delta(bt*Yk - deltaU + al*bt*Yk.^2 + dyy_yk ,...
        bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,secondDerivative) - dyy_yk;
    figure(4)
    mesh(x(9:end-8),y(9:end-8),residual(9:end-8,9:end-8)');
    title('Residual (no boundary values)')
    %mesh(x,y,z4')
    %mesh(x(points(1):points(2)),y(points(3):points(4)),z4')
    
    figure(15)
    mesh(x,y,residual');
    title('Residual (with all boundary values)')
    return;
    figure(5)
    plot(y,bigU(zeroX,:)/theta(end)) %(1+end)/2  
    xlabel('y')
    title('x==0 cross-section');
    figure(6) 
    plot(x,bigU(:,zeroY)/theta(end))   % ((x.^2)').*
    xlabel('x')
    title('y==0 cross-section');

    figure(7)
    plot(y,bigU(1,:)) %(1+end)/2
    xlabel('y')
    title('x=x_s_t cross-section');
    figure(8)
    plot(x,bigU(:,1))
    xlabel('x')
    title('y=y_s_t cross-section');

    figure(19);
    mesh(x,y,bigUTimeDer');
    xlabel('x');    ylabel('y');
    title('du/dt');
    xlabel('x');    ylabel('y');
    axis([x(1) x(end) y(1) y(end) -0.1 0.1]);
    colorbar;
    caxis([-0.001 0.001]);
    view(90,90);

    %figure(4);
    %mesh(x(ij:end),y(lo:end),dyuu');
    %xlabel('x');    ylabel('y');

    %================================
    %================================
    %*sqrt(1-c^2/bt)

    fig_ss11=figure(11);
    set(fig_ss11, 'OuterPosition', [0.0      	30.0        380.0     340.0]);
    %*sqrt(1-c^2)
    mesh(x,y,bigU');
    xlabel('x');    ylabel('y');
    title('end solution')
    %axis([x_st2 x_end2 y_st2 y_end2 -0.5 1]);
    axis([x_st x_end y_st y_end -0.01 .01]);
    colorbar;
    caxis([-0.00000001 .00000001]);
    view(0,90);

    fig_ss12=figure(12);
    set(fig_ss12, 'OuterPosition', [0.0      	30.0        380.0     340.0]);
    %*sqrt(1-c^2)
    mesh(x,y,( bigIC'));
    xlabel('x');    ylabel('y');
    title('IC')
    %axis([x_st2 x_end2 y_st2 y_end2 -0.5 1]);
    axis([x_st x_end y_st y_end -0.01 .01]);
    colorbar;
    caxis([-0.00000001 .00000001]);
    view(0,90);
end