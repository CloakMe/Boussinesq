clear;

for jl = 1:3
    yo(1,:) = '4';
    yo(2,:) = '2';
    yo(3,:) = '1';

    yo = cellstr(yo);
    ICType = 'Natali'; % Christov  Natali
    cellStr = strcat('ConvergenceTests\', ICType, 'IC_16_bt3_c045_h0', yo(jl), '_O(h^4)' );
    load (  cellStr{1} );
    %sum(tauVector)
    if(jl==1)
       vb = U(1:end,1:end); hb=h; 
    end
    if(jl==2)
       vc = U(1:end,1:end); hc=h;
    end
    if(jl==3)
       vd = U(1:end,1:end); 
    end
    PlotAssymptotics(x,y,h,zeroX,zeroY,bigU,0);
    clear('U');clear('cellStr');clear('yo');
end

clear('jl');
res_1 = vb - vc(1:2:end,1:2:end);
res_2 = vc - vd(1:2:end,1:2:end);

norm0402_L2 = hb*norm(res_1(:),2);
norm0201_L2 = hc*norm(res_2(:),2);
fprintf('||v_04 - v_02||_L2 = %.8e \n', norm0402_L2);
fprintf('||v_02 - v_01||_L2 = %.8e \n', norm0201_L2);
conv_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2);
fprintf('Conv_L2 = %.8e \n\n', conv_L2);

norm0402_Inf = max(max(abs(res_1(:))));
norm0201_Inf = max(max(abs(res_2(:))));
fprintf('||v_04 - v_02||_Inf = %.8e \n', norm0402_Inf);
fprintf('||v_02 - v_01||_Inf = %.8e \n', norm0201_Inf);
conv_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2);
fprintf('Conv_Inf = %.8e \n', conv_Inf);

clear;