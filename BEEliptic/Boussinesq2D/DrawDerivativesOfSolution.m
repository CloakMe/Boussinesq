function DrawDerivativesOfSolution(bigU,compBox,x,y,h,zeroX,zeroY,c,derivative,c1)

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

    sx = length(x)
    sy = length(y)
    mag = floor(20/h);
    xx=x(zeroX-mag+1:zeroX+mag-1); yy=y(zeroY-mag+1:zeroY+mag-1);
    %================================
    ox1 = ones(1,sy);
    X = x'*ox1;
    if(sx~=sy)
        oy1 = ones(1,sx);
        Y = (y'*oy1)';
     else
        Y = (y'*ox1)';
    end
    vdah = zeros(size(bigU));

    
    %================================

    figure(10);
    mag2 = floor(42);
    u_th = X.*YDer(bigU,derivative.first) - Y.*XDer(bigU,derivative.first) ;
    u_y = YDer(bigU,derivative.first) ;
    mesh(x,y,(u_th'));
    xlabel('x');    ylabel('y');
    title('u_t_h')
    axis([x_st x_end y_st y_end -0.1 .1]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.01 .01]);
    colorbar;
    caxis([-0.01 .01]);
    view(0,90);

    u_thth = Y.^2 .* XDer(bigU,derivative.second) +  X.^2 .* YDer(bigU,derivative.second) -...
        X.*XDer(bigU,derivative.first) - Y.*YDer(bigU,derivative.first) -...
        2*X.*Y.*XDer(YDer(bigU,derivative.first),derivative.first);
    figure(11);

    mesh(x*sqrt(1-c^2),y,(u_thth'));
    xlabel('x');    ylabel('y');
    title('u_t_h_t_h')
    axis([x_st x_end y_st y_end -0.001 0.001]);
    %axis([compBox.x_st2 x_end2 compBox.y_st2 compBox.y_end2 -0.001 0.001]);
    colorbar;
    caxis([-0.001 0.001]);
    view(0,90);






    vdah = zeros(size(bigU));

    u_yy =YDer(bigU,derivative.second);
    u_xx =XDer(bigU,derivative.second);
    u_yyyy =YDer(YDer(bigU,derivative.second),derivative.second);
    u_xxxx =XDer(XDer(bigU,derivative.second),derivative.second);
    u_xxyy =YDer(YDer(bigU,derivative.second),derivative.second);
    squareU_yy =YDer(bigU.^2,derivative.second);
    squareU_xx =XDer(bigU.^2,derivative.second);



    figure(1)
    mesh(x,y,(u_yy'));
    xlabel('x');    ylabel('y');
    title('u_y_y')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.0000001 0.0000001]);
    colorbar;
    caxis([-0.0000001 0.0000001]);
    view(90,90);

    figure(2)
    mesh(x,y,(u_xx'));
    xlabel('x');    ylabel('y');
    title('u_x_x')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.00000001 0.00000001]);
    colorbar;
    caxis([-0.00000001 0.00000001]);
    view(90,90);

    figure(3)
    mesh(x,y,(u_xxxx'));
    xlabel('x');    ylabel('y');
    title('u_x_x_x_x')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.000000001 0.000000001]);
    colorbar;
    caxis([-0.000000001 0.000000001]);
    view(90,90);

    figure(4)
    mesh(x,y,(u_yyyy'));
    xlabel('x');    ylabel('y');
    title('u_y_y_y_y')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.000000001 0.000000001]);
    colorbar;
    caxis([-0.000000001 0.000000001]);
    view(90,90);

    figure(5)
    mesh(x,y,(u_xxyy'));
    xlabel('x');    ylabel('y');
    title('u_x_x_y_y')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.0000000001 0.0000000001]);
    colorbar;
    caxis([-0.0000000001 0.0000000001]);
    view(90,90);

    figure(6)
    mesh(x,y,(squareU_xx'));
    xlabel('x');    ylabel('y');
    title('squ_x_x')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.00000000001 0.00000000001]);
    colorbar;
    caxis([-0.00000000001 0.00000000001]);
    view(90,90);

    figure(7)
    mesh(x,y,(squareU_yy'));
    xlabel('x');    ylabel('y');
    title('squ_y_y')
    axis([x_st x_end y_st y_end -0.000001 .000001]);
    %axis([compBox.x_st2 compBox.x_end2 compBox.y_st2 compBox.y_end2 -0.00000000001 0.00000000001]);
    colorbar;
    caxis([-0.00000000001 0.00000000001]);
    view(90,90);

    c12 = 1 - c^2;
    Boundary=c1*(c12*X.^2-Y.^2)./(c12*X.^2+Y.^2).^2;
    figure(8)
    mesh(x,y,Boundary' - bigU')
    axis([x_st x_end y_st y_end -0.0001 0.0001]);
    colorbar;
    caxis([-0.00001 0.00001]);
    view(0,90);