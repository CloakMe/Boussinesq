function [Qd,TF14] = pola2cart_v5(x, y, gama1,gama2, de_name, IC)
swit_ch = 1;
if(nargin == 6)
    swit_ch = 0;
end

[X,Y] = Domain(x,y);

lenY = length(y);
lenX = length(x);
   [zeroX,zeroY]=GetZeroNodes(x,y);
   if(zeroX==0 || zeroY == 0)
       error('Zero point not found! Suggest interval with zero point -> (0,0)');
   end
   if(zeroY>(lenY+1)/2 || (lenX+1)/2~=zeroX ) %|| ij>sy-lo+1
       error('Domain is disproportional! Check x-y boundaries!');
   end
   radiusEnd = floor(sqrt(x(end)^2 + y(end)^2)) + 1;
   if(swit_ch )       
       [r,dR,R,MAT]=nat_42(x,y,radiusEnd,10000,gama1,gama2, X, Y);
   else
       [r,dR,R,MAT]=BE_bkwrds(radiusEnd, gama1, de_name, IC, X, Y);
   end
       [valu_y ind_vy] = min(MAT);
       [valu indx] = min(valu_y);
       %MAT(indx,ind_vy(indx))
   
       TF = zeros(lenX,lenY);
       ii=-1;
       for i=zeroX:lenX
           ii=ii+1;
           for j = (zeroY+ii):lenY
             %[valu,ind] = max(MAT(i,j) == r);
             [valu,ind] = max(abs(MAT(i,j) - r) < 10^(-9) );
             TF(i ,j) = R(ind);
          end
       end
       
       TF14 = TF(zeroX:end,zeroY:end);
       if(lenY - zeroY >= lenX - zeroX)
           TF_mini = TF14(1:zeroX,1:zeroX);
           TF_mini = TF_mini + TF_mini' -diag(diag(TF_mini));
           TF14(1:end,1:lenX-zeroX+1)=TF_mini;
       else
           TF_mini = TF14(1:zeroY,1:zeroY);
           TF_mini = TF_mini + TF_mini' -diag(diag(TF_mini));
           TF14(1:lenY-zeroY+1,1:end)=TF_mini;
       end
       Qd = transf2qD(TF14,x,y,zeroX,zeroY);       
       
       %{
       T2 = TF + TF' -diag(diag(TF));
       T4 = T2((lenX+1)/2:lenX,(lenY+1)/2:lenY);
       T3 = T2;
       T3((1+end)/2,:) = 0;
       T2 = flipud(T2) + T3;
       T3 = T2;
       T3(:,(1+end)/2) = 0;
       TF=fliplr(T2) + T3;
       
       figure(1)
       mesh(x(zeroX:end),y(zeroY:end),TF(zeroX:end,zeroY:end)')
       figure(2)
       mesh(x,y,Qd')
       %}
end

function [X,Y]=Domain(x,y)
    
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
end
  
       
   
   