function PlotSolCharacteristics(btString, cString, hString ,orderString, domainLenString, solTypeString, bndString, legendString, additionalInfo) 
    
    len = length(legendString);
    if( length(btString) == 1 )
       btString = CreateMultiArrString( btString(1), len );
    end
    if( length(cString) == 1 )
       cString = CreateMultiArrString( cString(1), len );
    end
    if( length(hString) == 1 )
       hString = CreateMultiArrString( hString(1), len );
    end
    if( length(orderString) == 1 )
       orderString = CreateMultiArrString( orderString(1), len );
    end
    if( length(domainLenString) == 1 )
       domainLenString = CreateMultiArrString( domainLenString(1), len );
    end
    if( length(solTypeString) == 1 )
       solTypeString = CreateMultiArrString( solTypeString(1), len );
    end
    if( length(bndString) == 1 )
       bndString = CreateMultiArrString( bndString(1), len );
    end
    colors = { 'k', 'g', 'm', 'r', 'b', 'c'}; %
    figNum = 20;
    figure(figNum)
    for i=1:len
        bndCutSizeX = 0;
        bndCutSizeY = 0;
        shortBndString = '_ZB0';
        if( strcmp( bndString(i), 'ZeroBoundary' ) == true || strcmp( bndString(i), 'ZeroBnd' ) == true )
            shortBndString = '_ZB1';
        end
        cellStrEn = strcat('SavedWorkspaces\', bndString(i), '\', solTypeString(i), '_', domainLenString(i),...
            '_bt', btString(i), '_c0', cString(i),...
            '\', solTypeString(i), '_', domainLenString(i), shortBndString, '_bt',...
            btString(i), '_c0', cString(i), '_h0', hString(i), '_O(h^', orderString(i), ')' );

        warning('off','all');
        load (  cellStrEn{1} );
        warning('on','all');
        clear cellStrEn;
        
        current_hx = ( x(11) - x(1) ) / 10.0;
        current_hy = ( y(11) - y(1) ) / 10.0;
        bndPtsRemX = bndCutSizeX/current_hx;
        bndPtsRemY = bndCutSizeY/current_hy;
        
        figure(figNum)
        hold on;
        if( strcmp( additionalInfo, 'maximum' ) )
            ylabel('max|u_h|','FontSize',18);
            plot(t(1:end-1),max_v, colors{i} );    
            
        elseif( strcmp( additionalInfo, 'mass' ) )
            ylabel('Mass','FontSize',18);
            [indeces, shift] = BEUtilities.GetCommonIndexArray( t, II );
            indeces(1) = [];
            plot(t(indeces+shift),II(indeces),colors{i} );
            hold off;
            
        elseif( strcmp( additionalInfo, 'energy' ) )
            ylabel('E_h','FontSize',18);
            [indeces, shift] = BEUtilities.GetCommonIndexArray( t, EN );
            indeces(1) = [];
            plot( t(1:end-1),EN,colors{i} );            
        end
        hold off;
    end
    legend(legendString);    
    xlabel('t','FontSize',18);  
    set(gca,'FontSize',18);

end

function [multiArrString] = CreateMultiArrString( element, arrSize )
    multiArrString = { element };
    for i=1:arrSize
        multiArrString( i ) = element;
    end
end
