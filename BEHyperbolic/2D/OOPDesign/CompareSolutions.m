function CompareSolutions(btString, cString, hString ,orderString, additionalInfo) 

    fprintf('c = 0.%s, bt = %s, order = %s \n', cString, btString, orderString);
    
    [x1,y1,t1,EN1,II1,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString );
    [x2,y2,t2,EN2,II2,uEnTaylor] = GetBEEngineTaylorSol( btString, cString, hString, orderString );
    
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
        fprintf('tau enSave = %f\n', tau_enSave);
        tau_Taylor = (t2(end) - t2(1))/ (length( t2 ) - 1);
        fprintf('tau Taylor = %f\n', tau_Taylor);
    end
    
    if( additionalInfo == 2 )
        [minEN1, maxEN1, meanEN1] = GetValues(EN1);
        [minEN2, maxEN2, meanEN2] = GetValues(EN2);
        [minII1, maxII1, meanII1] = GetValues(II1);
        [minII2, maxII2, meanII2] = GetValues(II2);
        fprintf('         mean           min           max       \n');
        %fprintf('EnSave %4.6f %4.6f %4.6f\n', meanEN1, minEN1, maxEN1 );
        %fprintf('Taylor %4.6f %4.6f %4.6f\n', meanEN2, minEN2, maxEN2 );
        fprintf('EnSave %4.6f  & %4.6f  & %4.6f  \n', meanII1, minII1, maxII1 );
        fprintf('Taylor %4.6f  & %4.6f  & %4.6f  \n', meanII2, minII2, maxII2 );
        fprintf('----------------------------------------\n\n');
        return;
    end
    
    difference = (uEnTaylor - uEnSave);
    normDifference_L2 = h*norm(difference(:),2);
    fprintf('||v_Taylor - v_EnSave||_L2  = %.6f \n', normDifference_L2);

    normDifference_Inf = max(max(abs(difference(:))));
    fprintf('||v_Taylor - v_EnSave||_Inf = %.6f \n', normDifference_Inf);
    
    fprintf('||v_EnSave||_Inf            = %.6f \n', max(max(abs(uEnSave(:)))));
    fprintf('=========================\n\n');    
    
    if(additionalInfo == 0)
        return;
    end
    
    fprintf('sx = %d, sy = %d\n', length(x1), length(y1));    
    viewTypeX = 0;
    viewTypeY = 90;

    figure(14)
    mesh(x1,y1,uEnSave');
    view( viewTypeX, viewTypeY );
    colorbar;
    title('solution Energy Save');
    xlabel('x');            ylabel('y');

    figure(15)
    mesh(x2,y2,uEnTaylor');
    view( viewTypeX, viewTypeY );
    colorbar;
    title('solution Taylor');
    xlabel('x');            ylabel('y');

    figure(16)
    mesh(x1,y1,(uEnTaylor- uEnSave)');
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

