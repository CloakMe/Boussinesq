function [x,tt,max_v,EN,II,uEnSave] = GetBEEngineEnergySaveSol1D( tauSrting, hString, domainLen )
%v,dtv,va,tt,II
   
    if(nargin == 2)
        domainLen = '60';
        bndType = 0;
    end
    
    bndCutSizeX = 0;
        
    cellStrEn = strcat('SavedWorkspaces\Sol_vesi_O(tau^2 + h^2)_', domainLen, '_tau', tauSrting, '_h0', hString );

    warning('off','all');
    load (  cellStrEn );
    warning('on','all');

    current_hx = ( x(11) - x(1) ) / 10.0;
    bndPtsRemX = bndCutSizeX/current_hx;

    uEnSave = v( bndPtsRemX+1:end-bndPtsRemX );
    max_v = max(uEnSave);
    EN = 0;
end