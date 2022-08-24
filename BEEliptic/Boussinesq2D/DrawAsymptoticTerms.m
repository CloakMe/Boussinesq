function DrawAsymptoticTerms(x,y,h,al,bt,c,theta,bigU,bigUTimeDer,bigIC,U,compBox,secondDerivative, mu, extended)

    if(nargin == 14)
        extended = 0;
    end
    
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
    
    barF = floor(min(length(x), length(y))/20);
    
    [zeroX,zeroY]=GetZeroNodes(x,y);
    h = x(2) - x(1);
    lenx = floor(length( x(zeroX:end) )*2/5);
    leny = floor(length( y(zeroY:end) )*2/5);
    xpts = zeroX + lenx:length(x) - barF;
    ypts = zeroY + leny:length(y) - barF;
    xx = x(xpts);
    yy = y(ypts);
    Yk = bigU;
    domainUtils = BEDomainUtils( x, y, length(secondDerivative)-1, 2 ); 
    dyy_U = domainUtils.YDerivativeZeroBnd(Yk)/h^2;
    dxx_U = domainUtils.XDerivativeZeroBnd(Yk)/h^2;
    
    btU_xx = bt * dxx_U;
    bt1McSQU_yy = bt * (1-c^2) * dyy_U;
    
    coefU_yyyy = - (1 - bt * c^2) * domainUtils.YDerivativeZeroBnd( dyy_U )/h^2;
    coefU_xxyy = - (2 - bt * c^2) * domainUtils.YDerivativeZeroBnd( dxx_U )/h^2;
    coefU_xxxx =                  - domainUtils.XDerivativeZeroBnd( dxx_U )/h^2;
    
    coeffUsq_xx = + al * bt * domainUtils.XDerivativeZeroBnd(Yk.^2)/h^2;
    coeffUsq_yy = + al * bt * domainUtils.YDerivativeZeroBnd(Yk.^2)/h^2;
    
    residual = btU_xx + bt1McSQU_yy + coefU_xxxx + coefU_xxyy + coefU_yyyy + coeffUsq_xx + coeffUsq_yy;
    
    legendString = { '|\beta U_x_x|', '|\gamma_1U_y_y|', ...
        '|\gamma_2U_y_y_y_y|', '|\gamma_3U_x_x_y_y|', '|U_x_x_x_x|', ...
        '|\alpha\beta(U^2)_x_x|', '|\alpha\beta(U^2)_y_y|', '|R|', ...
         '|\gamma_4 y^{-4}|', '|\gamma_5 y^{-6}|'};
    
    figure(2)
    %ax2 = axes('Position',[0.1 0.1 0.64 0.80]);
    %plt = loglog(ax2, y(ypts), abs(btU_xx(zeroX, ypts)), 'b.',...
    plt = loglog(y(ypts), abs(btU_xx(zeroX, ypts)), 'b',...
         y(ypts), abs(bt1McSQU_yy(zeroX, ypts)), 'r:',...
         y(ypts), abs(coefU_yyyy(zeroX, ypts)), 'c', ...
         y(ypts), abs(coefU_xxyy(zeroX, ypts)), 'k',...
         y(ypts), abs(coefU_xxxx(zeroX, ypts)), ...
         y(ypts), abs(coeffUsq_xx(zeroX, ypts)),...
         y(ypts), abs(coeffUsq_yy(zeroX, ypts)), 'g',...
         y(ypts), abs(residual(zeroX, ypts) ), '--', ...
         y(ypts(1:5:end)), abs( ( btU_xx(zeroX, ypts(1))*y(ypts(1))^4 ) ./ y(ypts(1:5:end)).^4), 'ok', ...
         y(ypts(1:5:end)), abs( ( coefU_xxyy(zeroX, ypts(1))*y(ypts(1))^6 )./y(ypts(1:5:end)).^6), 'om');% y(zeroY+10:end-barF), 1 ./ y(zeroY+10:end-barF) .^ 2, 'k' ); %(1+end)/2  
     
    set(plt(5),'Color',[0.4500 0.120 0.90]);
    set(plt(6),'Color',[0.92 0.6 0.1]);    %orange
    set(plt(8),'Color',[0.20 0.550 0.20]);
    
    label_y = strcat('y in [', num2str(y(ypts(1)), '%.1f'), ',', num2str(y(ypts(end)), '%.1f'), ']');
    xlabel(label_y,'FontSize',17);    
    set(gca,'FontSize',18);
    legend(legendString, 'FontSize',14, 'Position',[0.772 0.166 0.227 0.829]);   %,'Location', 'eastoutside'
    title('x==0 cross-section','FontSize',16);
    %xlabel('x','FontSize',18);    ylabel('y','FontSize',18);
    
    
    legendString = {'|\gamma_6 x^{-4}|', '|\gamma_7 x^{-6}|'};
     
    figure(3) 
    plt = loglog(x(xpts(1:5:end)), abs( ( btU_xx(xpts(1), zeroY)*x(xpts(1))^4 ) ./ x(xpts(1:5:end)).^4), 'ok', ...
        x(xpts(1:5:end)), abs( ( coefU_xxyy(xpts(1), zeroY)*x(xpts(1))^6 )./x(xpts(1:5:end)).^6), 'om', ...
        x(xpts), abs(btU_xx(xpts, zeroY)), 'b',...
        x(xpts), abs(bt1McSQU_yy(xpts, zeroY)), 'r:',...
        x(xpts), abs(coefU_yyyy(xpts, zeroY)), 'c', ...
        x(xpts), abs(coefU_xxyy(xpts, zeroY)), 'k',...
        x(xpts), abs(coefU_xxxx(xpts, zeroY)), ...
        x(xpts), abs(coeffUsq_xx(xpts, zeroY)) ,...
        x(xpts), abs(coeffUsq_yy(xpts, zeroY)), 'g',...
        x(xpts), abs(residual(xpts, zeroY)),'--' );
    
    set(plt(5+2),'Color',[0.4500 0.120 0.90]);
    set(plt(6+2),'Color',[0.92 0.6 0.1]);    %255,158,0
    set(plt(8+2),'Color',[0.20 0.550 0.20]);
    xlim([x(xpts(1)) x(xpts(end))])
    label_x = strcat('x in [',num2str(x(xpts(1)),'%.1f'), ',  ', num2str(x(xpts(end)),'%.1f'), ']');
    xlabel(label_x,'FontSize',18);
    set(gca,'FontSize',18);
    legend(legendString, 'FontSize',14, 'Position',[0.205 0.2 0.15 0.2]);   %,'Location', 'eastoutside'
    title('y==0 cross-section','FontSize',16);    
    
    if(extended == 0)
        return;
    end
    
    figure(4)
    mesh(xx,yy, btU_xx(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('bt * U_x_x (no boundary values)');
    xlabel('x');    ylabel('y');    
    figure(5)
    mesh(xx,yy, bt1McSQU_yy(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('bt * (1-c^2) * U_y_y (no boundary values)');
    xlabel('x');    ylabel('y');    
    
    figure(10)
    mesh(xx,yy, coefU_yyyy(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title(' - (1 - bt * c^2) * U_y_y_y_y (no boundary values)');
    xlabel('x');    ylabel('y');
    figure(11)
    mesh(xx,yy, coefU_xxyy(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title(' - (2 - bt * c^2) * U_x_x_y_y (no boundary values)');
    xlabel('x');    ylabel('y');    
    figure(12)
    mesh(xx,yy, coefU_xxxx(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title(' - U_x_x_x_x (no boundary values)');
    xlabel('x');    ylabel('y');    
   
    figure(13)
    mesh(xx,yy, coeffUsq_xx(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('al * bt * (U^2)_x_x (no boundary values)')
    xlabel('x');    ylabel('y');        
    figure(14)
    mesh(xx,yy, coeffUsq_yy(zeroX + lenx:end-barF,zeroY - leny:zeroY + leny)');
    title('al * bt * (U^2)_y_y (no boundary values)')
    xlabel('x');    ylabel('y');
        
    figure(15)
    mesh(x,y, residual');
    title('Residual (no boundary values)');
    xlabel('x');    ylabel('y');
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