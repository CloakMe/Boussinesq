function DrawAsymptoticTerms(x,y,h,al,bt,c,theta,prmtrs,bigUTimeDer,bigIC,U,compBox,secondDerivative, mu, extended)

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
    
    barF = 0; %floor(min(length(x), length(y))/20);
    
    [zeroX,zeroY]=GetZeroNodes(x,y);

    lenx = floor(length( x(zeroX:end) )*2/5);
    leny = floor(length( y(zeroY:end) )*2/5);
    xpts = 1 + lenx:(length(x) + 1)/2 - barF;
    ypts = 1 + leny:(length(y) + 1)/2 - barF;
    xx = x(zeroX + lenx:end-barF);
    yy = y(zeroY + leny:end-barF);
    if(prmtrs.useZeroBoundary > 0 )
        domainUtils = BEDomainUtils( x, y, length(secondDerivative)-1, 2 ); 
        dyy_U = domainUtils.YDerivativeEvenFunZeroBnd(U)/h^2;
        dxx_U = domainUtils.XDerivativeEvenFunZeroBnd(U)/h^2;

        btU_xx = bt * dxx_U;
        bt1McSQU_yy = bt * (1-c^2) * dyy_U;

        coefU_yyyy = - (1 - bt * c^2) * domainUtils.YDerivativeEvenFunZeroBnd( dyy_U )/h^2;
        coefU_xxyy = - (2 - bt * c^2) * domainUtils.YDerivativeEvenFunZeroBnd( dxx_U )/h^2;
        coefU_xxxx =                  - domainUtils.XDerivativeEvenFunZeroBnd( dxx_U )/h^2;

        coeffUsq_yy = + al * bt * domainUtils.YDerivativeEvenFunZeroBnd(U.^2)/h^2;
        coeffUsq_xx = + al * bt * domainUtils.XerivativeEvenFunZeroBnd(U.^2)/h^2;
        
    else
        
        x2 = [x(zeroX:end) x(end)+h x(end)+2*h x(end)+3*h x(end)+4*h];
        y2 = [y(zeroY:end) y(end)+h y(end)+2*h y(end)+3*h y(end)+4*h];
        [X,Y]=Domain(x2,y2);
        
        c12 = 1-c^2; 
        approxBoundaryF = (c12 * X.^2 - Y.^2)./(c12*X.^2 + Y.^2).^2; 
        approxBoundaryF_xx =  6*c12* (Y.^4 - 6*c12*X.^2.*Y.^2 + c12 * X.^4) ./ (Y.^2 + c12*X.^2).^4;
        approxBoundaryF_yy = -6 *    (Y.^4 - 6*c12*X.^2.*Y.^2 + c12 * X.^4) ./ (Y.^2 + c12*X.^2).^4;
        
        outerTopBoundaryF_xx = approxBoundaryF_xx(1:end-4,end-3:end); 
        outerRigthBoundaryF_xx = approxBoundaryF_xx(end-3:end,1:end-4); 
        
        outerTopBoundaryF_yy = approxBoundaryF_yy(1:end-4,end-3:end); 
        outerRigthBoundaryF_yy = approxBoundaryF_yy(end-3:end,1:end-4); 
        
        outerTopBoundaryF = approxBoundaryF(1:end-4,end-3:end); 
        outerRigthBoundaryF = approxBoundaryF(end-3:end,1:end-4); 
        
        zeroMatrix = zeros(size(U));
        dyy_U = YDerivativeEvenFunctions(U, zeroMatrix, mu.muU*outerTopBoundaryF, secondDerivative);
        dxx_U = XDerivativeEvenFunctions(U, zeroMatrix, mu.muU*outerRigthBoundaryF, secondDerivative);

        btU_xx = bt * dxx_U;
        bt1McSQU_yy = bt * (1-c^2) * dyy_U;

        coefU_yyyy = - (1 - bt * c^2) * YDerivativeEvenFunctions(dyy_U, zeroMatrix, mu.muU*outerTopBoundaryF_yy, secondDerivative);
        coefU_xxyy = - (2 - bt * c^2) * YDerivativeEvenFunctions( dxx_U, zeroMatrix, mu.muU*outerTopBoundaryF_xx, secondDerivative );
        coefU_xxxx =                  - XDerivativeEvenFunctions( dxx_U, zeroMatrix, mu.muU*outerRigthBoundaryF_xx, secondDerivative );
        
        coeffUsq_yy = + al * bt * YDerivativeEvenFunctions(U.^2, zeroMatrix, mu.muU^2*outerTopBoundaryF.^2, secondDerivative);
        coeffUsq_xx = + al * bt * XDerivativeEvenFunctions(U.^2, zeroMatrix, mu.muU^2*outerRigthBoundaryF.^2, secondDerivative);
    end
    
    residual = btU_xx + bt1McSQU_yy + coefU_xxxx + coefU_xxyy + coefU_yyyy + coeffUsq_xx + coeffUsq_yy;
    
    legendString = { '|\beta U_x_x|', '|\gamma_1U_y_y|', ...
        '|\gamma_2U_y_y_y_y|', '|\gamma_3U_x_x_y_y|', '|U_x_x_x_x|', ...
        '|\alpha\beta(U^2)_x_x|', '|\alpha\beta(U^2)_y_y|', '|R|', ...
         '|\gamma_4 y^{-4}|', '|\gamma_5 y^{-6}|'};
    
     %( btU_xx(zeroX, ypts(1))*yy(1)^4 )
     %
    figure(2)
    %ax2 = axes('Position',[0.1 0.1 0.64 0.80]);
    %plt = loglog(ax2, yy, abs(btU_xx(zeroX, ypts)), 'b.',...
    plt = loglog(yy, abs(btU_xx(1, ypts)), 'b',...
         yy, abs(bt1McSQU_yy(1, ypts)), 'r:',...
         yy, abs(coefU_yyyy(1, ypts)), 'c', ...
         yy, abs(coefU_xxyy(1, ypts)), 'k',...
         yy, abs(coefU_xxxx(1, ypts)), ...
         yy, abs(coeffUsq_xx(1, ypts)),...
         yy, abs(coeffUsq_yy(1, ypts)), 'g',...
         yy, abs(residual(1, ypts) ), '--', ...
         yy(1:5:end), abs( 6 * bt * (1-c^2) * mu.muU ./ yy(1:5:end).^4), 'ok', ...
         yy(1:5:end), abs( ( coefU_xxyy(1, ypts(1))*yy(1)^6 ) ./ yy(1:5:end).^6), 'om');% y(zeroY+10:end-barF), 1 ./ y(zeroY+10:end-barF) .^ 2, 'k' ); %(1+end)/2  
     
    set(plt(5),'Color',[0.4500 0.120 0.90]);
    set(plt(6),'Color',[0.92 0.6 0.1]);    %orange
    set(plt(8),'Color',[0.20 0.550 0.20]);
    
    label_y = strcat('y in [', num2str(y(ypts(1)), '%.1f'), ',', num2str(yy(end), '%.1f'), ']');
    xlabel(label_y,'FontSize',17);    
    set(gca,'FontSize',18);
    legend(legendString, 'FontSize',14, 'Position',[0.772 0.166 0.227 0.829]);   %,'Location', 'eastoutside'
    title('x==0 cross-section','FontSize',16);
    %xlabel('x','FontSize',18);    ylabel('y','FontSize',18);
    
    
    legendString = {'|\gamma_6 x^{-4}|', '|\gamma_7 x^{-6}|'};
     %( btU_xx(xpts(1), zeroY)*xx(1)^4 )
     %
    figure(3) 
    plt = loglog(xx(1:5:end), abs( 6 * bt * mu.muU / (1-c^2) ./ xx(1:5:end).^4), 'ok', ...
        xx(1:5:end), abs( ( coefU_xxyy(xpts(1), 1)*xx(1)^6 ) ./xx(1:5:end).^6), 'om', ...
        xx, abs(btU_xx(xpts, 1)), 'b',...
        xx, abs(bt1McSQU_yy(xpts, 1)), 'r:',...
        xx, abs(coefU_yyyy(xpts, 1)), 'c', ...
        xx, abs(coefU_xxyy(xpts, 1)), 'k',...
        xx, abs(coefU_xxxx(xpts, 1)), ...
        xx, abs(coeffUsq_xx(xpts, 1)) ,...
        xx, abs(coeffUsq_yy(xpts, 1)), 'g',...
        xx, abs(residual(xpts, 1)),'--' );
    
    set(plt(5+2),'Color',[0.4500 0.120 0.90]);
    set(plt(6+2),'Color',[0.92 0.6 0.1]);    %255,158,0
    set(plt(8+2),'Color',[0.20 0.550 0.20]);
    xlim([xx(1) xx(end)])
    label_x = strcat('x in [',num2str(xx(1),'%.1f'), ',  ', num2str(xx(end),'%.1f'), ']');
    xlabel(label_x,'FontSize',18);
    set(gca,'FontSize',18);
    legend(legendString, 'FontSize',14, 'Position',[0.205 0.2 0.15 0.2]);   %,'Location', 'eastoutside'
    title('y==0 cross-section','FontSize',16);    
    
    if(extended == 0)
        return;
    end
    
    figure(4)
    mesh(xx,yy, btU_xx(1 + lenx:end-barF,1:1 + leny)');
    title('bt * U_x_x (no boundary values)');
    xlabel('x');    ylabel('y');    
    figure(5)
    mesh(xx,yy, bt1McSQU_yy(1 + lenx:end-barF,1:1 + leny)');
    title('bt * (1-c^2) * U_y_y (no boundary values)');
    xlabel('x');    ylabel('y');    
    
    figure(10)
    mesh(xx,yy, coefU_yyyy(1 + lenx:end-barF,1:1 + leny)');
    title(' - (1 - bt * c^2) * U_y_y_y_y (no boundary values)');
    xlabel('x');    ylabel('y');
    figure(11)
    mesh(xx,yy, coefU_xxyy(1 + lenx:end-barF,1:1 + leny)');
    title(' - (2 - bt * c^2) * U_x_x_y_y (no boundary values)');
    xlabel('x');    ylabel('y');    
    figure(12)
    mesh(xx,yy, coefU_xxxx(1 + lenx:end-barF,1:1 + leny)');
    title(' - U_x_x_x_x (no boundary values)');
    xlabel('x');    ylabel('y');    
   
    figure(13)
    mesh(xx,yy, coeffUsq_xx(1 + lenx:end-barF,1:1 + leny)');
    title('al * bt * (U^2)_x_x (no boundary values)')
    xlabel('x');    ylabel('y');        
    figure(14)
    mesh(xx,yy, coeffUsq_yy(1 + lenx:end-barF,1:1 + leny)');
    title('al * bt * (U^2)_y_y (no boundary values)')
    xlabel('x');    ylabel('y');
        
    figure(15)
    mesh(xx,yy, residual');
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