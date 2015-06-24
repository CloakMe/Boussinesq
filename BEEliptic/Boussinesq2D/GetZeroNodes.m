function [zeroX,zeroY]=GetZeroNodes(x,y)
    sx = length(x);
    sy = length(y);
    zeroX = 0;zeroY = 0;
    for k = 1:sx
       if( -10^(-11)<x(k) && x(k)<10^(-11))
           zeroX=k;break;
       end
   end
   for k = 1:sy
       if( -10^(-11)<y(k) && y(k)<10^(-11))
           zeroY=k;break;
       end
   end
   if(zeroX == 0 || zeroY == 0)
       error('Cannot find Mid Point!');
   end
end