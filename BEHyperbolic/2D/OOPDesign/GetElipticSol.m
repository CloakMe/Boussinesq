function [bigUcut, c] = GetElipticSol( btString, cString, hString, orderString, withBoundary, domainLen )

    if(nargin == 5)
        domainLen = '40';
    end
    
    bndCutSizeX = 0;
    bndCutSizeY = 0;
    bndString = 'WithBoundary';
    if( withBoundary == 0 )
        bndString = 'ZeroBoundary';
    end
    cellStrTlr = strcat('..\..\..\BEEliptic\Boussinesq2D\', bndString, '\ChristovIC_', domainLen, '_bt', btString, '_c0', cString,...
            '\Oh', orderString, '\ChristovIC_', domainLen, '_ZB1_bt',...
            btString, '_c0', cString, '_h0', hString, '_O(h^', orderString, ')' );
        
    warning('off','all');
    load (  cellStrTlr );
    warning('on','all');
    
    current_hx = ( x(end) - x(1) ) / ( size(x,1) * size(x,2) );
    current_hy = ( y(end) - y(1) ) / ( size(y,1) * size(y,2) );
    bndPtsRemX = bndCutSizeX/current_hx;
    bndPtsRemY = bndCutSizeY/current_hy;
    
    bigUcut = bigU( bndPtsRemX+1:end-bndPtsRemX, bndPtsRemY+1:end-bndPtsRemY );
end