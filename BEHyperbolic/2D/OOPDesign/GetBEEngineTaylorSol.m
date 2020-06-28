function [x,y,t,max_v,EN,II,uEnTaylor] = GetBEEngineTaylorSol( btString, cString, hString, orderString, withBoundary, domainLen )

    if(nargin == 5)
        domainLen = '40';
    end
    
    bndCutSizeX = 0;
    bndCutSizeY = 0;
    
    if( withBoundary == 0 )
        cellStrTlr = strcat('SavedWorkspaces\ZeroBoundary\Taylor_', domainLen, '_bt', btString, '_c0', cString,...
            '\Taylor_', domainLen, '_ZB1_bt',...
            btString, '_c0', cString, '_h0', hString, '_O(h^', orderString, ')' );
    elseif ( withBoundary == 1 )
        cellStrTlr = strcat('SavedWorkspaces\WithBoundary\Taylor_', domainLen, '_bt', btString, '_c0', cString,...
            '\Taylor_', domainLen, '_ZB0_bt',...
            btString, '_c0', cString, '_h0', hString, '_O(h^', orderString, ')' );
    end
    warning('off','all');
    load (  cellStrTlr );
    warning('on','all');
    
    current_hx = ( x(end) - x(1) ) / ( size(x,1) * size(x,2) );
    current_hy = ( y(end) - y(1) ) / ( size(y,1) * size(y,2) );
    bndPtsRemX = bndCutSizeX/current_hx;
    bndPtsRemY = bndCutSizeY/current_hy;
    
    uEnTaylor = vl( bndPtsRemX+1:end-bndPtsRemX, bndPtsRemY+1:end-bndPtsRemY );
end