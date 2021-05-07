function CompareTaylorVsEnergySave(btString, cString, hString ,orderString, domainLen, additionalInfo) 

    fprintf('c = 0.%s, bt = %s, order = %s \n', cString, btString, orderString);
    
    [x1,y1,t1,max_v1,EN1,II1,uEnSave] = GetBEEngineEnergySaveSol( btString, cString, hString, domainLen );
    [x2,y2,t2,max_v2,EN2,II2,uEnTaylorZeroBoundary] = GetBEEngineTaylorSol( btString, cString, hString, orderString, 0, domainLen );
    
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
    
    if( additionalInfo == 1 || additionalInfo == 2 )
        difference = (uEnTaylorZeroBoundary - uEnSave);
        normDifference_L2 = h*norm(difference(:),2);
        fprintf('||v_Taylor - v_EnSave||_L2  = %.3e \n', normDifference_L2);

        normDifference_Inf = max(max(abs(difference(:))));
        fprintf('||v_Taylor - v_EnSave||_Inf = %.3e \n', normDifference_Inf);

        fprintf('||v_EnSave||_Inf            = %.6f \n', max(max(abs(uEnSave(:)))));
        fprintf('=========================\n\n');    

        if(additionalInfo == 1)
            return;
        end
        figNumber = 1;
        for i = 1:50
            figNumber = i;
            if(  ishandle(i) == false )
                break;
            end
        end
        viewTypeX = 0;
        viewTypeY = 90;
        figure(figNumber)
        [zeroX,zeroY]=GetZeroNodes(x1,y1);
        magX = floor( (x1(end)-x1(zeroX)) / (3*h) );
        magY = floor( (y1(end)-y1(zeroY)) / (3*h) );
        xIndeces = zeroX-magX+1:zeroX+magX-1;
        yIndeces = zeroY-magY+1:zeroY+magY-1;
        xx=x1(xIndeces); yy=y1(yIndeces); 
        mesh(xx,yy,(uEnSave(xIndeces,yIndeces) - uEnTaylorZeroBoundary(xIndeces,yIndeces))');
        view( viewTypeX, viewTypeY );
        colorbar;
        title('solution difference');
        xlabel('x');            ylabel('y');

        figure(figNumber+1)
        plot(t1(2:end),max_v1, 'r', t2(2:end),max_v2, 'b');
        title('solution maximums');
        xlabel('t');            %ylabel('solution maximums');
        legend('EnergySave','Taylor');
    elseif( additionalInfo == 3 || additionalInfo == 4 )

        [minEN1, maxEN1, meanEN1] = GetValues(EN1);
        [minEN2, maxEN2, meanEN2] = GetValues(EN2);
        [minII1, maxII1, meanII1] = GetValues(II1);
        [minII2, maxII2, meanII2] = GetValues(II2);
        fprintf('========================================\n');  
        fprintf('Energy   min          max        |diff|       \n');
        fprintf('EnSave %4.6f & %4.6f & %4.6f\n', minEN1, maxEN1, maxEN1 - minEN1 );
        fprintf('Taylor %4.6f & %4.6f & %4.6f\n', minEN2, maxEN2, maxEN2 - minEN2 );
        fprintf('----------------------------------------\n');
        fprintf('Integral   min        max        |diff|       \n');
        fprintf('EnSave: %4.6f  & %4.6f  & %4.6f  \n', minII1, maxII1, maxII1 - minII1 );
        fprintf('Taylor: %4.6f  & %4.6f  & %4.6f  \n', minII2, maxII2, maxII2 - minII2 );
        fprintf('----------------------------------------\n\n');
    
        if(additionalInfo == 3)
            return;
        end
        
        figNumber = 1;
        for i = 1:50
            figNumber = i;
            if(  ishandle(i) == false )
                break;
            end
        end
        figure(figNumber)        
        %hold on;
        st1 = round( length(t1)/length(II1) );   
        newt1 = t1(st1:st1:end);
        [indeces1, shift1] = BEUtilities.GetCommonIndexArray( newt1, II1 );
        indeces1(1) = [];
        
        st2 = round( length(t2)/length(II2) );
        newt2 = t2(st2:st2:end);
        [indeces2, shift2] = BEUtilities.GetCommonIndexArray( newt2, II2 );
        indeces2(1) = [];
        
        plot(newt1(indeces1+shift1),II1(indeces1),'r',newt1(1),II1(2)+II1(2)/1000.0,newt1(end),II1(end)-II1(end)/1000.0 , ...
            newt2(indeces2+shift2),II2(indeces2),'b',newt2(1),II2(2)+II2(2)/1000.0,newt2(end),II2(end)-II2(end)/1000.0 )
        title('Integral EnergySave/Taylor');
        legend('EnergySave','Taylor');
        xlabel('time "t"');  ylabel('I');
        %hold off;
    end
  
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

