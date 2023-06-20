function CompareStartEndSolitons(btString, cString, hString ,orderString, domainLen, additionalInfo, color, tFix) 

    if( nargin == 6 )
        color = 'k';
    end
	if( nargin == 7 )
        tFix = false;
    end
%c = 0.9
%P:\PhDWork\OneSoliton_bt1\ChristovIC_40_80_bt1_c090_h010_O(h^6).mat
%P:\PhDWork\OneSoliton_bt1\Hyperb_40_80_bt1_c090_h010_O(h^6).mat
%clear
%c = 0.52;
%load "C:\Boussinesq\BEHyperbolic\2D\OOPDesign\SavedWorkspaces\Hyperb_40_bt3_c052_h005_O(h^6).mat"
%load "C:\Boussinesq\BEEliptic\Boussinesq2D\WithBoundary\ChristovIC_40_bt3_c052\Oh6\ChristovIC_40_bt3_c052_h005_O(h^6).mat"

    
    
    %get solution from elliptic problem.
    [bigU, c] = GetElipticSol( btString, cString, hString, orderString, 1, domainLen );
    %get solution from hyperbolic problem
    [x,y,t,max_v,EN,II,vl] = GetBEEngineTaylorSol( btString, cString, hString, orderString, 1, domainLen, additionalInfo, tFix );
    
    fprintf('c = 0.%s, bt = %s, order = %s, h = %f, tau = %f \n', cString, btString, orderString, x(2)-x(1), t(2)-t(1));
    
    if( additionalInfo == 2 )
        stepM = 1;%ceil(length(max_v)/40);
        mySlice = 1:stepM:length(max_v);
        %if(orderString == '4')
        %    stepM = length(max_v)/10;
        %    mySlice = 1:stepM:length(max_v);
        %end
        if(mySlice(end) ~= length(max_v) && stepM > 1)
            mySlice = [mySlice length(max_v)];
        end
        figure(17)
        hold on;
        plot(t(mySlice), max_v(mySlice), color ); %t(1:stepM:end-1)
        hold off;
        %title('Evolution of the maximum');
        xlabel('t','FontSize',16);  ylabel('max(v_h)','FontSize',16); 
        set(gca,'FontSize',16);
        return;
    end
    

    
    viewTypeX = 90;
    viewTypeY = 90;

    x_idx = GetIdx( x, 0 );
    y_start_idx = GetIdx( y, 0 );
    
    y_end_idx = GetIdx( y, t(end)*c );
    
    y_half_size = min( length(y) - y_end_idx - 20, floor( (y_end_idx - y_start_idx)/2 ) );
    x_half_size = length(x) - x_idx - 20;

    indeces_x = x_idx - x_half_size:x_idx + x_half_size;
    start_indeces_y = y_start_idx - y_half_size:y_start_idx + y_half_size;
    end_indeces_y = y_end_idx - y_half_size:y_end_idx + y_half_size;

    %return;
    x_vec = x(indeces_x);
    y_start = y(start_indeces_y);
    y_end = y(end_indeces_y);

    vl_start = bigU( indeces_x, start_indeces_y );
    vl_end = vl( indeces_x, end_indeces_y );

%     figure(1)
%     mesh(x_vec, y_start, vl_start');
%     title(['Solution at time: ', num2str(0)]);
%     xlabel('x'); ylabel('y');
%     colorbar;
%     axis tight;
%     view( viewTypeX, viewTypeY );
% 
%     figure(2)
%     mesh(x_vec, y_end, vl_end');
%     title(['Solution at time: ', num2str(t(end))]);
%     xlabel('x'); ylabel('y');
%     colorbar;
%     axis tight;
%     view( viewTypeX, viewTypeY );
% 
%     figure(3)
%     mesh(x_vec, y_start, (vl_start-vl_end)');
%     title('Solution difference between start and end time');
%     xlabel('x'); ylabel('y');
%     colorbar;
%     axis tight;
%     view( viewTypeX, viewTypeY );
% 
%     figure(4)
%     mesh(x, y, vl');
%     title('Solution');
%     xlabel('x'); ylabel('y');
%     colorbar;
%     axis tight;    
    
    vl_diff = (vl_start-vl_end);
    vl_diff_L2 = (x(2) - x(1))*norm( vl_diff(:), 2 );
    %fprintf('||vl_diff||_L2 = %.6f \n', vl_diff_L2);
    
    vl_diff_Inf = max(max(abs(vl_diff(:))));
    fprintf('||vl_diff||_L2, ||vl_diff||_Inf = %.6f & %.6f \n', vl_diff_L2, vl_diff_Inf);
    

end