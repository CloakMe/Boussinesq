function [x,y,t,max_v,EN,II,uEnTaylor] = GetBEEngineTaylorSol( btString, cString, hString, orderString, withBoundary, domainLen, additionalInfo )

    if(nargin == 5)
        domainLen = '40';
    end
    
   if(nargin <= 6)
		additionalInfo = 1;
    end

    bndCutSizeX = 0;
    bndCutSizeY = 0;
    domainLenInName = domainLen;
    if( additionalInfo ~= 2 && ...
        ( strcmp(domainLen, '128') == 1 && strcmp(btString, '1') == 1 && strcmp(hString, '40') == 1 && strcmp(cString, '90') == 1) ...
            || (strcmp(domainLen, '30') == 1 && strcmp(btString, '3') == 1 && strcmp(hString, '20') == 1 && strcmp(cString, '45') == 1) )
        %domainLenInName = strcat(domainLen, '_TFix');
    end
    if(strcmp(hString, '40') == 1 && additionalInfo ~= 2)
        newHstring = strcat(hString, '_Tfix');
    else
        newHstring = hString;
    end
    if( withBoundary == 0 )
        cellStrTlr = strcat('SavedWorkspaces\ZeroBoundary\Taylor_', domainLen, '_bt', btString, '_c0', cString,...
            '\Taylor_', domainLenInName, '_ZB1_bt',...
            btString, '_c0', cString, '_h0', newHstring, '_O(h^', orderString, ')' );
    elseif ( withBoundary == 1 )
        cellStrTlr = strcat('SavedWorkspaces\WithBoundary\Taylor_', domainLen, '_bt', btString, '_c0', cString,...
            '\Taylor_', domainLenInName, '_ZB0_bt',...
            btString, '_c0', cString, '_h0', newHstring, '_O(h^', orderString, ')' );
    elseif ( withBoundary == 2 )
        cellStrTlr = strcat('SavedWorkspaces\ZeroBnd\TaylorZeroBnd_', domainLen, '_bt', btString, '_c0', cString,...
            '\TaylorZeroBnd_', domainLenInName, '_ZB2_bt',...
            btString, '_c0', cString, '_h0', newHstring, '_O(h^', orderString, ')' );
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