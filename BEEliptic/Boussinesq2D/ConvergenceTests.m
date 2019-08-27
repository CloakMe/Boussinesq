clear;
remPnts = 0;
for jl = 1:3
    yo(1,:) = '20';
    yo(2,:) = '10';
    yo(3,:) = '05';

    yo = cellstr(yo);
    ICType = 'Christov'; % Christov  Natali
    cellStr = strcat('SavedWorkspaces\', ICType, 'IC_40_bt1_c090_h0', yo(jl), '_O(h^2)' );
    load (  cellStr{1} );
    %sum(tauVector)
    if(jl==1)
       vb = U(1+remPnts:end - remPnts,1+remPnts:end - remPnts); hb=h; 
       size(U)
    end
    if(jl==2)
        remPnts = 2*remPnts;
       vc = U(1+remPnts:end - remPnts,1+remPnts:end - remPnts); hc=h;
       size(U)
    end
    if(jl==3)
        remPnts = 2*remPnts;
       vd = U(1+remPnts:end - remPnts,1+remPnts:end - remPnts); 
       size(U)
    end
    %PlotAssymptotics(x,y,h,zeroX,zeroY,bigU,0);
    clear('U');clear('cellStr');clear('yo');
end

clear('jl');
res_1 = vb - vc(1:2:end,1:2:end);
res_2 = vc - vd(1:2:end,1:2:end);

norm0402_L2 = hb*norm(res_1(:),2);
norm0201_L2 = hc*norm(res_2(:),2);
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