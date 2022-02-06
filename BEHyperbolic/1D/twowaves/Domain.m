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