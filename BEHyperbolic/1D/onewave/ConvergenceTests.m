clear;
bndCutSizeX = 0;
bndCutSizeY = 0;

bndDirectory = ''; % 'ZeroBoundary\' 'WithBoundary\' ZeroBnd
nameN = 'Sol_Taylor_v3_NoMixedDer'; % 'Sol_Taylor_v3_' 'Taylor_' TaylorZeroBnd_
domLength = '100_'; %'80_' '40_', '128_'
orderString = '4';

%=======================================================
% 
% yo(1,:) = 'tau0000500_h020';
% yo(2,:) = 'tau0000250_h020';
% yo(3,:) = 'tau0000125_h020';

yo(1,:) = 'h020';
yo(2,:) = 'h010';
yo(3,:) = 'h005';

yo = cellstr(yo);

for jl = 1:3

    cellStr = strcat('SavedWorkspaces\', nameN, '_O(tau^', orderString, ' + h^', orderString, ')_', domLength, yo(jl) );
    %cellStr = strcat('Hyp_40_bt1_c090_h0', yo(jl), '_O(h^4)' );
    warning('off','all');
    load (  cellStr{1} );
    warning('on','all');
    fprintf('tau = %f, h = %f \n', tau, h);

    bndPtsRemX = 0.0;
    bndPtsRemY = 0.0;
    
    if(jl==1)
       vb = v; 
       hb=h;
       taub = tau;
    end
    if(jl==2)
       vc = v;
       hc=h;
       tauc = tau;
    end
    if(jl==3)
       vd = v;
       hd=h;
       taud = tau;
    end
    clear('v');clear('cellStr');
end

% E1 = norm(v_h1 - v_h05(1:2:end),2)
% E2 = norm(v_h05(1:2:end) - v_h025(1:4:end),2)
% 
% conv = (log(E1)-log(E2))/log(2)

clear('jl');
if(hb == hc && hd == hb)
    res_1 = (vb - vc);
    res_2 = (vc - vd);
else
    res_1 = (vb - vc(1:2:end));
    res_2 = (vc - vd(1:2:end));
end

norm0402_L2 = hb*norm(res_1(:),2);
norm0201_L2 = hc*norm(res_2(:),2);
fprintf('Solution Convergence:\n');

fprintf('||v_04 - v_02||_L2 = %.6f \n', norm0402_L2);
fprintf('||v_02 - v_01||_L2 = %.6f \n', norm0201_L2); %%.8e
conv_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2);
fprintf('Conv_L2 = %.8e \n\n', conv_L2);

norm0402_Inf = max(max(abs(res_1(:))));
norm0201_Inf = max(max(abs(res_2(:))));
fprintf('||v_04 - v_02||_Inf = %.6f \n', norm0402_Inf);
fprintf('||v_02 - v_01||_Inf = %.6f \n', norm0201_Inf);
conv_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2);
fprintf('Conv_Inf = %.8e \n', conv_Inf);



clear;