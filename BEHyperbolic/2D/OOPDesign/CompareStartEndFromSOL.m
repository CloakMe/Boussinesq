function CompareStartEndFromSOL(h, c, step, endStep, SOLPath, vl_start) 

%c = 0.9
%P:\PhDWork\OneSoliton_bt1\ChristovIC_40_80_bt1_c090_h010_O(h^6).mat
%P:\PhDWork\OneSoliton_bt1\Hyperb_40_80_bt1_c090_h010_O(h^6).mat
%clear
%c = 0.52;
%load "C:\Boussinesq\BEHyperbolic\2D\OOPDesign\SavedWorkspaces\Hyperb_40_bt3_c052_h005_O(h^6).mat"
%load "C:\Boussinesq\BEEliptic\Boussinesq2D\WithBoundary\ChristovIC_40_bt3_c052\Oh6\ChristovIC_40_bt3_c052_h005_O(h^6).mat"

    fprintf('h = %f \n', h);
    %get solution from hyperbolic problem
    for iter=step:step:endStep
        
        cellStrTlr = strcat(SOLPath, '\vz_',  int2str(iter), '.mat' );
        warning('off','all');
        load (  cellStrTlr );
        vl_end = vz;
        warning('on','all');
        
        distance = iter*c;
        shiftSteps = distance / h;
        vl_end_smaller = vl_end(:,shiftSteps+1:end);
        
        vl_start_smaller = vl_start(:,1:end-shiftSteps);
        vl_diff = (vl_start_smaller-vl_end_smaller);
        vl_diff_L2 = h*norm( vl_diff(:), 2 );
        %fprintf('||vl_diff||_L2 = %.6f \n', vl_diff_L2);
        
        vl_diff_Inf = max(max(abs(vl_diff(:))));
        fprintf('iter = %d \n', iter);
        fprintf('||vl_diff||_L2, ||vl_diff||_Inf = %.6f & %.6f \n', vl_diff_L2, vl_diff_Inf);
    end


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
    

    

end