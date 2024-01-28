clear;
bndCutSizeX = 0;
bndCutSizeY = 0;

bndDirectory = ''; % 'ZeroBoundary\' 'WithBoundary\' ZeroBnd
nameN = 'Sol_Taylor_v3_NoMixedDer'; % 'Sol_Taylor_v3_' 'Taylor_' TaylorZeroBnd_
domLength = '100_'; %'80_' '40_', '128_'
tauOrderString = '4';
hOrderString = '4';
%=======================================================
% 
% yo(1,:) = 'tau0000500_h020';
% yo(2,:) = 'tau0000250_h020';
% yo(3,:) = 'tau0000125_h020';

hString = ''; %'h020'
if(isempty(hString))
    yo(1,:) = 'h020';
    yo(2,:) = 'h010';
    yo(3,:) = 'h010';
else
    yo(1,:) = 'tau0001000_';
    yo(2,:) = 'tau0000500_';
    yo(3,:) = 'tau0000500_';
end

tauString = 'tau0000010_';

yo = cellstr(yo);
useExactSolution = 1;

endTime = 0;
x_stN = 0;
x_endN = 0;
for jl = 1:3

    try
        if(isempty(hString))
            cellStr = strcat('SavedWorkspaces\', nameN, '_O(tau^', tauOrderString, ' + h^', hOrderString, ')_', domLength, tauString, yo(jl) );
        else
            cellStr = strcat('SavedWorkspaces\', nameN, '_O(tau^', tauOrderString, ' + h^', hOrderString, ')_', domLength, yo(jl), hString ); 
        end
        warning('off','all');
        load (  cellStr{1} );
        warning('on','all');
        fprintf('tau = %f, h = %f \n', tau, h);
    catch me
        fprintf('Could not load requested file %s!\n', cellStr{1});
    end
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
       x_stN = start_x; 
       x_endN = end_x;
       endTime = t_start + t_interval;
    end
    if(jl==3)
        if(~isempty(hString))
            taud = tauc/2
            hd = h;
        else
            hd=hc/2
            taud = tau;
        end
        if(useExactSolution == 1)
            ic_utils = IC_2Waves();
            x = x_stN:hd:x_endN;
            vd = ic_utils.GetInitialCondition2w(x, endTime)';
        else
            vd = v; 
        end
    end
    clear('v');clear('cellStr');
end

% E1 = norm(v_h1 - v_h05(1:2:end),2)
% E2 = norm(v_h05(1:2:end) - v_h025(1:4:end),2)
% 
% conv = (log(E1)-log(E2))/log(2)

clear('jl');
if(isempty(hString))
    res_1 = (vb - vc(1:2:end));
    res_2 = (vc - vd(1:2:end));
else
    res_1 = (vb - vc);
    res_2 = (vc - vd);
end

if(~isempty(hString))
    norm0402_L2 = norm(res_1(:),2);
    norm0201_L2 = norm(res_2(:),2);
else
    norm0402_L2 = hb*norm(res_1(:),2);
    norm0201_L2 = hc*norm(res_2(:),2);    
end
fprintf('Solution Convergence:\n');

fprintf('||v_04 - v_02||_L2 = %.10f \n', norm0402_L2);
fprintf('||v_02 - v_01||_L2 = %.10f \n', norm0201_L2); %%.8e
conv_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2);
fprintf('Conv_L2 = %.8e \n\n', conv_L2);

norm0402_Inf = max(max(abs(res_1(:))));
norm0201_Inf = max(max(abs(res_2(:))));
fprintf('||v_04 - v_02||_Inf = %.10f \n', norm0402_Inf);
fprintf('||v_02 - v_01||_Inf = %.10f \n', norm0201_Inf);
conv_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2);
fprintf('Conv_Inf = %.8e \n', conv_Inf);



clear;