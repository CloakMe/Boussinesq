function [x,y,t,EN,II,uEnTaylor] = GetBEEngineTaylorSol( btString, cString, hString, orderString, bndFlag )

    bndCutSizeX = 0;
    bndCutSizeY = 0;
    
    if( bndFlag == 0 )
        cellStrTlr = strcat('SavedWorkspaces\ZeroBoundary\Taylor_40_ZB1_bt',...
            btString, '_c0', cString, '_h0', hString, '_O(h^', orderString, ')' );
    elseif ( bndFlag == 1 )
        cellStrTlr = strcat('SavedWorkspaces\WithBoundary\Hyperb_40_ZB0_bt',...
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