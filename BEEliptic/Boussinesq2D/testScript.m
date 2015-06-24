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
    
    mag = floor(5/h);
    xx=x(zeroX-mag+1:zeroX+mag-1); yy=y(zeroY-mag+1:zeroY+mag-1); 

    BigUNearCenter = bigU(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);
    bigUTimeDerNearCenter = bigUTimeDer(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);
    % bigICNearCenter = bigIC(zeroX-mag+1:zeroX+mag-1,zeroY-mag+1:zeroY+mag-1);


    figure(2);
    mesh(xx,yy,( bigUTimeDerNearCenter'))%bbU'
    xlabel('x');    ylabel('y');
    title('time derivative')
    xx1 = 0:0.4:16;
    xx2 = 0:0.2:16;
    xx3 = 0:0.1:16;
    figure(5);
    mesh(xx1,xx1,(res_1)');
    xlabel('x');    ylabel('y');
    title('difference')
    
    figure(6);
    mesh(xx2,xx2,(res_2)');
    xlabel('x');    ylabel('y');
    title('difference')
    
    
    norm0402_L2 = hb*norm(res_1(:),2);
    norm0201_L2 = hc*norm(res_2(:),2);
    fprintf('||v_04 - v_02||_L2 = %.8e \n', norm0402_L2);
    fprintf('||v_02 - v_01||_L2 = %.8e \n', norm0201_L2);
    conv_L2 = log(abs(norm0201_L2/norm0402_L2))/log(2);
    fprintf('Conv_L2 = %.8e \n\n', conv_L2);