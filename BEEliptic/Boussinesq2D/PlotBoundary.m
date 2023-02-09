function PlotBoundary(x, y, U, approxBoundaryF, domainUtils, al, bt, c, residual)
    
    %[yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(length(x),length(y),step);
    
    outerTopBoundaryF=approxBoundaryF(1:end-4,end-3:end); 
    outerRigthBoundaryF=approxBoundaryF(end-3:end,1:end-4); 
    h = x(2)-x(1);
    dxx_U = domainUtils.XDerivativeEvenFunZeroBnd(U)/h^2;
    dyy_U = domainUtils.YDerivativeEvenFunZeroBnd(U)/h^2;
        
    btU_xx = bt * dxx_U;
    bt1McSQU_yy = bt * (1-c^2) * dyy_U;
    
    coefU_yyyy = - (1 - bt * c^2) * domainUtils.YDerivativeEvenFunZeroBnd( dyy_U )/h^2;
    coefU_xxyy = - (2 - bt * c^2) * domainUtils.YDerivativeEvenFunZeroBnd( dxx_U )/h^2;
    coefU_xxxx =                  - domainUtils.XDerivativeEvenFunZeroBnd( dxx_U )/h^2;
    
    coeffUsq_xx = + al * bt * domainUtils.XDerivativeEvenFunZeroBnd(U.^2)/h^2;
    coeffUsq_yy = + al * bt * domainUtils.YDerivativeEvenFunZeroBnd(U.^2)/h^2;
    
    [zeroX,zeroY]=GetZeroNodes(x,y);
    %lenx = floor(length( x(zeroX:end) )*2/5);
    leny = floor(length( y(zeroY:end) )*2/5);
    ypts = 1 + leny:(length(y) + 1)/2;
    %hy = y(2)-y(1);
    figure(4)
    plt = loglog(y(ypts), abs(btU_xx(zeroX, ypts)), 'b',...
         y(ypts), abs(bt1McSQU_yy(zeroX, ypts)), 'r:',...
         y(ypts), abs(coefU_yyyy(zeroX, ypts)), 'm+', ...
         y(ypts), abs(coefU_xxyy(zeroX, ypts)), 'bo', ...
         y(ypts), abs(coefU_xxxx(zeroX, ypts)), 'r', ...
         y(ypts), abs(coeffUsq_xx(zeroX, ypts)), 'yo', ...
         y(ypts), abs(coeffUsq_yy(zeroX, ypts)), 'g',...
         y(ypts), abs(residual(zeroX, ypts) ), '--'); %, ...
         %y(ypts(1:5:end)), abs( 6 * mu.muU ./ y(ypts(1:5:end)).^4), 'ok', ...
         %y(ypts(1:5:end)), abs( ( coefU_xxyy(zeroX, ypts(1))*y(ypts(1))^6 ) ./ y(ypts(1:5:end)).^6), 'om');% y(zeroY+10:end-barF), 1 ./ y(zeroY+10:end-barF) .^ 2, 'k' ); %(1+end)/2  
     
     legendString = { '|\beta U_x_x|', '|\gamma_1U_y_y|', ...
        '|\gamma_2U_y_y_y_y|', '|\gamma_3U_x_x_y_y|', '|U_x_x_x_x|', ...
        '|\alpha\beta(U^2)_x_x|', '|\alpha\beta(U^2)_y_y|', '|R|'};%, ...
         %'|\gamma_4 y^{-4}|', '|\gamma_5 y^{-6}|'};
     legend(legendString, 'FontSize',14, 'Position',[0.772 0.166 0.227 0.829]);   %,'Location', 'eastoutside'
%     figure(4);
%     mesh(x(xLen), y(yHeiA), abs( U(xLen, yHeiA) - muU*approxBoundaryF(xLen, yHeiA) )' );
%     xlabel('x');    ylabel('y');
%     view(0,90);
%     colorbar;
%     
%     figure(5);
%     mesh(x(xHeiA), y(yLen), abs( U(xHeiA, yLen) - muU*approxBoundaryF(xHeiA, yLen) )' );
%     xlabel('x');    ylabel('y');
%     view(0,90);
%     colorbar;
end
    %x2 = [x x(end)+h x(end)+2*h x(end)+3*h x(end)+4*h];
    %y2 = [y y(end)+h y(end)+2*h y(end)+3*h y(end)+4*h];
    
    %mesh(x(zeroX:end),y(zeroX:end),U');