clear;
ICType = 'Christov'; % Christov  Natali

for jl = 1:3
    for jp=1:2
        yo(1,:) = '4';
        yo(2,:) = '2';
        yo(3,:) = '1';

        yo = cellstr(yo);
        
        if(jp==1)
            ICType = 'Christov';
        else
            ICType = 'Natali';
        end
        cellStr = strcat('ConvergenceTests\', ICType, 'IC_16_bt3_c045_h0', yo(jl), '_O(h^4)' );
        load (  cellStr{1} );
        %sum(tauVector)
        if(jl==1)
           if(jp==1)
               vb = U(1:end,1:end); hb=h; 
           else
               vb = vb - U(1:end,1:end);
               PlotAssymptotics(x,y,h,zeroX,zeroY,vb,1);
           end
        end
        if(jl==2)
           if(jp==1)
               vc = U(1:end,1:end); hc=h;
           else
               vc = vc - U(1:end,1:end); 
               PlotAssymptotics(x,y,h,zeroX,zeroY,vc,1);
           end
        end
        if(jl==3)
           if(jp==1)
               vd = U(1:end,1:end); hd=h;
           else
               vd = vd - U(1:end,1:end);
               PlotAssymptotics(x,y,h,zeroX,zeroY,vd,1);
           end           
        end
        clear('U');clear('cellStr'); clear('yo');
    end
    
   
end

clear('jl');
clear('jp');

res_b = vb;
res_c = vc;
res_d = vd;

normhb_L2 = hb*norm(res_b(:),2);
normhc_L2 = hc*norm(res_c(:),2);
normhd_L2 = hd*norm(res_d(:),2);
fprintf('||v_Chr - v_Nat||_b_L2 = %.8e \n', normhb_L2);
fprintf('||v_Chr - v_Nat||_c_L2 = %.8e \n', normhc_L2);
fprintf('||v_Chr - v_Nat||_d_L2 = %.8e \n\n', normhd_L2);

convbc_L2 = log(abs(normhb_L2/normhc_L2))/log(2);
convcd_L2 = log(abs(normhc_L2/normhd_L2))/log(2);
fprintf('Conv_bc_L2 = %.8e \n', convbc_L2);
fprintf('Conv_cd_L2 = %.8e \n\n', convcd_L2);

normhb_Inf = max(max(abs(res_b(:))));
normhc_Inf = max(max(abs(res_c(:))));
normhd_Inf = max(max(abs(res_d(:))));

fprintf('||v_Chr - v_Nat||_b_L2 = %.8e \n', normhb_Inf);
fprintf('||v_Chr - v_Nat||_c_L2 = %.8e \n', normhc_Inf);
fprintf('||v_Chr - v_Nat||_d_L2 = %.8e \n\n', normhd_Inf);

convbc_Inf = log(abs(normhb_Inf/normhc_Inf))/log(2);
convcd_Inf = log(abs(normhc_Inf/normhd_Inf))/log(2);
fprintf('conv_bc_Inf = %.8e \n', convbc_Inf);
fprintf('conv_cd_Inf = %.8e \n', convcd_Inf);
clear;