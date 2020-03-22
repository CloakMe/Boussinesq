function [x,y,t,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString )

    bndCutSizeX = 0;
    bndCutSizeY = 0;

    cellStrEn = strcat('SavedWorkspaces\ZeroBoundary\EnergySave_40_ZB1_bt',...
        btString, '_c0', cString, '_h0', hString, '_O(h^2)' );

    warning('off','all');
    load (  cellStrEn );
    warning('on','all');

    current_hx = ( x(11) - x(1) ) / 10.0;
    current_hy = ( y(11) - y(1) ) / 10.0;
    bndPtsRemX = bndCutSizeX/current_hx;
    bndPtsRemY = bndCutSizeY/current_hy;

    uEnSave = vl( bndPtsRemX+1:end-bndPtsRemX, bndPtsRemY+1:end-bndPtsRemY );
    
end