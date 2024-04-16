function [x,tt,max_v,EN,II,uEnTaylor] = GetBEEngineTaylorSol1D( tauSrting, hString, domainLen )

    if(nargin == 2)
        domainLen = '60';
    end
      
    bndCutSizeX = 0;
    
    cellStrTlr = strcat('SavedWorkspaces\Sol_Taylor_v3_org_O(tau^2 + h^2)_', domainLen, '_tau', tauSrting, '_h0', hString );
    warning('off','all');
    load (  cellStrTlr );
    warning('on','all');
    
    current_hx = ( x(end) - x(1) ) / ( size(x,1) * size(x,2) );
    bndPtsRemX = bndCutSizeX/current_hx;
    
    uEnTaylor = v( bndPtsRemX+1:end-bndPtsRemX );
    max_v = max(uEnTaylor);
    EN = 0;
end