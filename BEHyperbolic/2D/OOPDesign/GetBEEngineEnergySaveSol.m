function [x,y,t,max_v,EN,II,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString, domainLen, bndType )

    if(nargin == 4)
        bndType = 0;
    end
    
    if(nargin == 3)
        domainLen = '40';
        bndType = 0;
    end
    
    if(bndType ==0)
        bndDir = 'WithBoundary';
        bndTypeString = '_ZB0';
    elseif(bndType == 1)
        bndDir = 'ZeroBoundary';
        bndTypeString = '_ZB1';
    else
        bndDir = 'ZeroBnd';
        bndTypeString = '_ZB2';
    end
    bndCutSizeX = 0;
    bndCutSizeY = 0;
        
    cellStrEn = strcat('SavedWorkspaces\', bndDir, '\EnergySave_', domainLen, '_bt', btString, '_c0', cString,...
        '\EnergySave_', domainLen, bndTypeString, '_bt',...
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