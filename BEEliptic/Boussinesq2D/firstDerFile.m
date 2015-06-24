function firstDerFile()
    for gg=[4 2 1]
        clearvars -except gg;
        %TestConvIterMethod  ||  SavedSolutions
        load (['TestConvIterMethod\who_30_bt3_c045_h0' num2str(gg) '_ord6']);
        firstDer=GetFiniteDifferenceCoeff([0,1,2],1)/h;
        U_y = U(1,1:3)*firstDer;
        fprintf('U_y = %.4e \n', U_y);
        
        gBIC = generateIC(x,y,c,bt1);
        IC_y = gBIC(zeroX,zeroY:zeroY+2)*firstDer;
        fprintf('IC_y = %.4e \n', IC_y);
    end
end

function [gBIC] = generateIC(x,y,c,bt1)
  sx = length(x);
  sy = length(y);
  % HA4AJLHu YcJLOBuQ _____________
   ox1 = ones(1,sy);
   X = x'*ox1;
   if(sx~=sy)
        oy1 = ones(1,sx);
        Y = (y'*oy1)';
   else
        Y = (y'*ox1)';
   end
   
   for k = 1:sx
       if( -10^(-11)<x(k) && x(k)<10^(-11))
           ij=k;break;
       end
   end
   for k = 1:sy
       if( -10^(-11)<y(k) && y(k)<10^(-11))
           lo=k;break;
       end
   end
   gBIC = u_ex2d_mat_vc(X,Y,c,bt1);%/sqrt(bt1)
end