function CompareTaylorWithBoundaryVsZeroBoundary(btString, cString, hString ,orderString, additionalInfo) 

    fprintf('c = 0.%s, bt = %s, order = %s \n', cString, btString, orderString);
    
    %[x1,y1,t1,EN1,II1,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString );
    [x1,y1,t1,EN1,II1,uEnTaylorWithBoundary] = GetBEEngineTaylorSol( btString, cString, hString, orderString, 1 );
    [x2,y2,t2,EN2,II2,uEnTaylorZeroBoundary] = GetBEEngineTaylorSol( btString, cString, hString, orderString, 0 );
    
    if( length( x1 ) ~= length( x2 ) || length( y1 ) ~= length( y2 ) )
        fprintf('Different sizes in X, Y or T - dimensions!\n');
        return;
    elseif( sum( x1 ~= x2 ) && sum( y1 ~= y2 ) && t1(end) ~= t2(end) )
        fprintf('Different values in x, y, t dimension vectors!\n');
        return;
    end

    h = (x1(end) - x1(1))/ (length( x1 ) - 1);
    tau_enSave = (t1(end) - t1(1))/ (length( t1 ) - 1);
    fprintf('h   = %f\n',  h);
    if( length( t1 ) == length( t2 ) )
        fprintf('tau = %f\n', tau_enSave); 
    else
        fprintf('tau Taylor With Boundary = %f\n', tau_enSave);
        tau_Taylor = (t2(end) - t2(1))/ (length( t2 ) - 1);
        fprintf('tau Taylor Zero Boundary = %f\n', tau_Taylor);
    end
    
    if( additionalInfo == 2 )
        [minEN1, maxEN1, meanEN1] = GetValues(EN1);
        [minEN2, maxEN2, meanEN2] = GetValues(EN2);
        [minII1, maxII1, meanII1] = GetValues(II1);
        [minII2, maxII2, meanII2] = GetValues(II2);
        fprintf('         mean           min           max       \n');
        %fprintf('EnSave %4.6f %4.6f %4.6f\n', meanEN1, minEN1, maxEN1 );
        %fprintf('Taylor %4.6f %4.6f %4.6f\n', meanEN2, minEN2, maxEN2 );
        fprintf('Taylor WB%4.6f  & %4.6f  & %4.6f  \n', meanII1, minII1, maxII1 );
        fprintf('Taylor ZB%4.6f  & %4.6f  & %4.6f  \n', meanII2, minII2, maxII2 );
        fprintf('----------------------------------------\n\n');
        return;
    end
    
    difference = (uEnTaylorZeroBoundary - uEnTaylorWithBoundary);
    normDifference_L2 = h*norm(difference(:),2);
    %fprintf('||v_Taylor - v_EnSave||_L2  = %.6f \n', normDifference_L2);
    fprintf('||v_Taylor_WB - v_Taylor_ZB||_L2  = %.6f \n', normDifference_L2);

    normDifference_Inf = max(max(abs(difference(:))));
    %fprintf('||v_Taylor - v_EnSave||_Inf = %.6f \n', normDifference_Inf);
    fprintf('||v_Taylor_WB - v_Taylor_ZB||_Inf = %.6f \n', normDifference_Inf);
    
    %fprintf('||v_EnSave||_Inf            = %.6f \n', max(max(abs(uEnSave(:)))));
    fprintf('||v_Taylor_WB||_Inf            = %.6f \n', max(max(abs(uEnTaylorWithBoundary(:)))));
    fprintf('=========================\n\n');    
    
    if(additionalInfo == 0)
        return;
    end
        
    fprintf('sx = %d, sy = %d\n', length(x1), length(y1));    
    viewTypeX = 0;
    viewTypeY = 90;

    if(additionalInfo == 3)
        figure(18)        
        %hold on;
        [indeces, shift] = BEUtilities.GetCommonIndexArray( t1, II1 );
        indeces(1) = [];
        %plot(t1(indeces+shift),II1(indeces),'k',t1(1),II1(2)+II1(2)/1000.0,t1(end),II1(end)-II1(end)/1000.0 )
        %
        st1 = round( length(t1)/length(II1) );
        plot(t1(st1:st1:end),II1,'k',t1(1),II1(2)+II1(2)/1000.0,t1(end),II1(end)-II1(end)/1000.0 );
        %hold off;
        title('Integral Taylor With Boundary');
        xlabel('time "t"');  ylabel('I');
        
        figure(19)
        %hold on;
        [indeces, shift] = BEUtilities.GetCommonIndexArray( t2, II2 );
        indeces(1) = [];
        plot(t2(indeces+shift),II2(indeces),'k',t2(1),II2(2)+II2(2)/1000.0,t2(end),II2(end)-II2(end)/1000.0 )
        %hold off;
        title('Integral Taylor Zero Boundary');
        xlabel('time "t"');  ylabel('I');
        return;
    end
    
    figure(14)
    mesh(x1,y1,uEnTaylorWithBoundary');
    view( viewTypeX, viewTypeY );
    colorbar;
    title('solution Taylor With Boundary');
    xlabel('x');            ylabel('y');

    figure(15)
    mesh(x2,y2,uEnTaylorZeroBoundary');
    view( viewTypeX, viewTypeY );
    colorbar;
    title('solution Taylor Zero Boundary');
    xlabel('x');            ylabel('y');

    figure(16)
    mesh(x1,y1,(uEnTaylorZeroBoundary- uEnTaylorWithBoundary)');
    view( viewTypeX, viewTypeY );
    colorbar;
    title('solution difference');
    xlabel('x');            ylabel('y');
end

function [minValue, maxValue, meanValue] = GetValues(array)
    if( length( array ) <= 1 )
       breakHere = 1; 
    end
    if( array( 1 ) == 0 )
        array = array(2:end);
    end
    minValue = min(array);
    maxValue = max(array);
    meanValue = sum(array)/length(array);
end

