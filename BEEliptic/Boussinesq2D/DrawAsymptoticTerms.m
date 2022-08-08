function DrawAsymptoticTerms(x,y,h,al,bt,c,theta,bigU,bigUTimeDer,bigIC,U,compBox,secondDerivative)
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
    
    barF = 100;
    [zeroX,zeroY]=GetZeroNodes(x,y);
    lenx = floor(length( x(zeroX:end) )*1/2);
    leny = floor(length( y ) / 4);
    xx = x(zeroX + lenx:end-barF);
    yy = y(zeroY - leny:zeroY + leny);
    Yk = bigU;
    dyy_U = YDer(Yk,secondDerivative);
    dxx_U = XDer(Yk,secondDerivative);
    res0 = bt * ((1-c^2) * dyy_U + dxx_U);
    res0Y = bt * ((1-c^2) * dyy_U);
    res0X = bt * (dxx_U);
    
    figure(4)
    mesh(xx,yy, res0(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('bt * ( (1-c^2) * U_y_y + U_x_x) (no boundary values)');
    
    figure(41)
    mesh(xx,yy, res0Y(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('bt * ((1-c^2) * U_y_y ) (no boundary values)');
    
    figure(42)
    mesh(xx,yy, res0X(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('bt * ( U_x_x) (no boundary values)');
    
%     figure(15)
%     plot(y(zeroY:end-barF), bigU(zeroX, zeroY:end-barF), 'b', y(zeroY:end-barF), res(zeroX, zeroY:end-barF), 'r', y(zeroY+10:end-barF), 1 ./ y(zeroY+10:end-barF) .^ 2, 'k' ); %(1+end)/2  
%     xlabel('y')
%     title('x==0 cross-section');
%     figure(16) 
%     plot(x(zeroX:end-barF), bigU(zeroX:end-barF, zeroY), 'b', x(zeroX:end-barF), res(zeroX:end-barF, zeroY), 'r', x(zeroX+10:end-barF), 1 ./ x(zeroX+10:end-barF) .^ 2, 'k')   % ((x.^2)').*
%     xlabel('x')
%     title('y==0 cross-section');
    YkSquare = Yk.^2;
    res1 = + al * bt * ( YDer(YkSquare,secondDerivative) + XDer(YkSquare,secondDerivative) );
    res1Y = + al * bt * ( YDer(YkSquare,secondDerivative) );
    
    figure(5)
    mesh(xx,yy, res1(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('al * bt * Delta(U^2) (no boundary values)')
    
    figure(51)
    mesh(xx,yy, res1Y(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('al * bt * (U^2)_y_y (no boundary values)')
    
    YkSquare = Yk.^2;
    res1X = + al * bt * ( XDer(YkSquare,secondDerivative) );
    figure(52)
    mesh(xx,yy, res1X(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('al * bt * (U^2)_x_x (no boundary values)')
    
    res2 = YDer( (bt * c^2 - 1) * dyy_U - dxx_U, secondDerivative ) + XDer( (bt * c^2 - 1) * dyy_U - dxx_U, secondDerivative );
    figure(6)
    mesh(xx,yy, res2(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('Delta( (bt * c^2 - 1)U_y_y - U_x_x )(no boundary values)');

    figure(11)
    mesh(xx,yy, (res0(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny) + res1(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny) +...
        + res2(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny))');
    title('Residual (no boundary values)');
    
return;
    bigZeroMatrix=zeros(size(bigU));
    dyy_yk = bt*c^2*YDer(Yk,secondDerivative);
    deltaU = Delta(Yk,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,secondDerivative);
    residual =  Delta(bt*Yk - deltaU + al*bt*Yk.^2 + dyy_yk ,...
        bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,secondDerivative) - dyy_yk;
    
    ff0 = ( - dyy_yk + bt * YDer(Yk, secondDerivative) + bt * XDer(Yk, secondDerivative) ) ;
    %res1 = + al * bt * ( YDer(Yk.^2,secondDerivative) + XDer(Yk.^2,secondDerivative) );
    ff1 = Delta( al*bt*Yk.^2,...
        bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,secondDerivative);
    
    %res2 = YDer( (bt * c^2 - 1) * dyy_U - dxx_U, secondDerivative ) + XDer( (bt * c^2 - 1) * dyy_U + dxx_U, secondDerivative );
    
%     ff21 = YDer(- deltaU + dyy_yk, secondDerivative);
%     ff22 = XDer(- deltaU + dyy_yk, secondDerivative);
    ff2 = Delta(- deltaU + dyy_yk ,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,bigZeroMatrix,secondDerivative);
    figure(45)
    mesh(x,y,(ff0+ff1+ff2)');
    title('Residual (with all boundary values)')
    
    figure(55)
    mesh(x,y,(ff0-res0)');
    figure(56)
    mesh(x,y,(ff1-res1)');
    figure(57)
    mesh(x,y,(ff2-res2)');
end