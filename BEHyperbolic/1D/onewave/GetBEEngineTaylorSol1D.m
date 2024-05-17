function [x,tt,max_v,EN,II,uEnTaylor, t_start, t_interval] = GetBEEngineTaylorSol1D( tauSrting, hString, domainLen, orderT )

    if(nargin == 3)
        orderT = 2;
    end
    if(nargin == 2)
        domainLen = '60';
        orderT = 2;
    end
      
    bndCutSizeX = 0;
    cellStrTlr = strcat('SavedWorkspaces\Sol_Taylor_v3_org_O(tau^2 + h^2)_', domainLen, '_tau', tauSrting, '_h0', hString );
    if(orderT == 4)
        cellStrTlr = strcat('SavedWorkspaces\Sol_Taylor_v3_org_O(tau^4 + h^4)_', domainLen, '_tau', tauSrting, '_h0', hString );
    end
    warning('off','all');
    load (  cellStrTlr );
    warning('on','all');
    
    current_hx = ( x(end) - x(1) ) / ( size(x,1) * size(x,2) );
    bndPtsRemX = bndCutSizeX/current_hx;
    
    uEnTaylor = v( bndPtsRemX+1:end-bndPtsRemX );
    max_v = max(va);
    EN = 0;
end