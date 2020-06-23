function [x,y,t,max_v,EN,II,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString, domainLen )

    if(nargin == 3)
        domainLen = '40';
    end
    
    bndCutSizeX = 0;
    bndCutSizeY = 0;
        
    cellStrEn = strcat('SavedWorkspaces\ZeroBoundary\EnergySave_', domainLen, '_bt', btString, '_c0', cString,...
        '\EnergySave_', domainLen, '_ZB1_bt',...
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