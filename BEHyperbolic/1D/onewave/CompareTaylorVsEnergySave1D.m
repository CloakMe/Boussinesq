function CompareTaylorVsEnergySave1D(tauString, hString, domainLen, additionalInfo, orderT) 

    if(nargin == 4)
        orderT = 2;
    end
    fprintf('tau = 0.%s, h = %s, domLen = %s \n', tauString, hString, domainLen);
    [x1,t1,max_v1,EN1,II1,uEnSave, t_start1, t_interval1] = GetBEEngineEnergySaveSol1D( tauString, hString, domainLen );
    [x2,t2,max_v2,EN2,II2,uTaylorZB, t_start2, t_interval2] = GetBEEngineTaylorSol1D( tauString, hString, domainLen, orderT );
   
    if( length( x1 ) ~= length( x2 ) )
        fprintf('Different sizes in X or T - dimensions!\n');
        return;
    elseif( (sum( x1 ~= x2 ) && t1(end) ~= t2(end)) )
        fprintf('Different values in x, y, t dimension vectors!\n');
        return;
    elseif( sum( t_start1 ~= t_start2 ) || t_interval1 ~= t_interval2 )
        fprintf('t_start1 = %f, t_interval1 = %f\n', t_start1, t_interval1);
        fprintf('t_start2 = %f, t_interval2 = %f\n', t_start2, t_interval2);
        fprintf('Different start times or time intervals!\n');
        return;        
    end
    
    ic_utils = IC_2Waves();
    u_end = ic_utils.GetInitialCondition2w(x1, t_start1+t_interval1)';
    
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
    figure(11)
    plot(x1, uTaylorZB, 'g', x1(1:5:end), u_end(1:5:end), 'b');
    title('End solution');
    if( additionalInfo == 1 || additionalInfo == 2 )
        
        difference = (u_end - uTaylorZB);
        normDifference_L2 = h*norm(difference(:),2);
        normDifference_Inf = max(max(abs(difference(:))));
        fprintf('|u_{ex}-u_Taylor|_L2 = %6e\n', normDifference_L2);
        fprintf('|u_{ex}-u_Taylor|_Inf = %6e\n', normDifference_Inf);        difference = (u_end - uEnSave);
        normDifference_L2 = h*norm(difference(:),2);
        normDifference_Inf = max(max(abs(difference(:))));
        fprintf('|u_{ex}-u_EnSave|_L2 = %6e\n', normDifference_L2);
        fprintf('|u_{ex}-u_EnSave|_Inf = %6e\n', normDifference_Inf);

        fprintf('||u_{ex}||_Inf            = %.6e \n', max(max(abs(u_end(:)))));
   
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
        zeroX=(length(x1)+1)/2;
        magX = floor( (x1(end)-x1(zeroX)) / (3*h) );
        xIndeces = zeroX-magX+1:zeroX+magX-1;
        xx=x1(xIndeces); 
        plot(xx,(uEnSave(xIndeces) - uTaylorZB(xIndeces))');
        xlabel('x','FontSize',18);   
        set(gca,'FontSize',18);

        figure(figNumber+1)
        plot(t1,max_v1, 'r', t2,max_v2, 'b');
        title('solution maximums');
        xlabel('t');            %ylabel('solution maximums');
        legend('EnergySave','Taylor');
        

    elseif( additionalInfo == 3 || additionalInfo == 4 )

        %[minEN1, maxEN1, meanEN1] = GetValues(EN1);
        %[minEN2, maxEN2, meanEN2] = GetValues(EN2);
        [minII1, maxII1, meanII1] = GetValues(II1);
        [minII2, maxII2, meanII2] = GetValues(II2);
        %fprintf('========================================\n');  
        %fprintf('Energy   min          max        |diff|       \n');
        %fprintf('EnSave %4.6f & %4.6f & %4.6f\n', minEN1, maxEN1, maxEN1 - minEN1 );
        %fprintf('Taylor %4.6f & %4.6f & %4.6f\n', minEN2, maxEN2, maxEN2 - minEN2 );
        fprintf('----------------------------------------\n');
        fprintf('Mass       min        max        |diff|       \n');
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
           
        plot(newt1(indeces1+shift1),II1(indeces1),'r^',newt2(indeces2+shift2),II2(indeces2),'bv', ...
            newt1(1),II1(2)+II1(2)/1000.0,newt1(end),II1(end)-II1(end)/1000.0, newt2(1),II2(2)+II2(2)/1000.0,newt2(end),II2(end)-II2(end)/1000.0 )
        title('Mass');
        legend('Conservative Scheme','Taylor');
        xlabel('time "t"');  ylabel('Mass');
        
        %
        if( ( h == 0.05 && btString == '3' ) || ( h == 0.1 && btString == '1' ) )
            colorPalette = 'b';
            if( orderString == '4' )
                colorPalette = 'g';
            elseif( orderString == '6' )
                colorPalette = 'r';
            end
            figure(30)
            hold on;       

            plot(newt2(indeces2+shift2), II2(indeces2),colorPalette )
            %, ...
            %    newt2(1),II2(2)+II2(2)/1000.0,newt2(end),II2(end)-II2(end)/1000.0
            title('Mass, Taylor');
            xlabel('time "t"');  ylabel('Mass');
            if( orderString == '6' )
                legend('O(h^2 + {\tau}^2)', 'O(h^4 + {\tau}^4)', 'O(h^6 + {\tau}^6)');
            end
            hold off;    
        end
        
        %
        [indeces1, shift1] = BEUtilities.GetCommonIndexArray( newt1, EN1 );
        indeces1(1) = [];
        
        [indeces2, shift2] = BEUtilities.GetCommonIndexArray( newt2, EN2 );
        indeces2(1) = [];
        
        figure(figNumber+10)
        plot(newt1(indeces1+shift1),EN1(indeces1),'r^',newt2(indeces2+shift2),EN2(indeces2),'bv', ...
            newt1(1),EN1(2)+EN1(2)/1000.0,newt1(end),EN1(end)-EN1(end)/1000.0, newt2(1),EN2(2)+EN2(2)/1000.0,newt2(end),EN2(end)-EN2(end)/1000.0 )
        
        title('Energy');  
        legend('Conservative Scheme','Taylor');
        xlabel('time "t"');  ylabel('Energy');
        %hold off;
        %
        if( ( h == 0.05 && btString == '3' ) || ( h == 0.1 && btString == '1' ) )
            colorPalette = 'b';
            if( orderString == '4' )
                colorPalette = 'g';
            elseif( orderString == '6' )
                colorPalette = 'r';
            end
            figure(31)
            hold on;       

            plot(newt2(indeces2+shift2),EN2(indeces2),colorPalette )
            %, ...
            %    newt2(1),II2(2)+II2(2)/1000.0,newt2(end),II2(end)-II2(end)/1000.0
            title('Energy, Taylor');
            xlabel('time "t"');  ylabel('Energy');
            if( orderString == '6' )
                legend('O(h^2 + {\tau}^2', 'O(h^4 + {\tau}^4', 'O(h^6 + {\tau}^6');
            end
            hold off;
        end
        
        %
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

