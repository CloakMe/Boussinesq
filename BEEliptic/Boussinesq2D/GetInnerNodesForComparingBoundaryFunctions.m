function [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step)
    
% The approximation boundary function (derived from the simplified BPE)
% and the current solution (with respect
% to the current iteration) are to be compared on the following node numbers:
% (xHeiA,yLen) and (xLen,yHeiA)
   band = floor(min(xLength, yLength)/4);
   
   stY = mod(yLength , step);
   if(stY == 0)
       stY = step;
   end
   yLen = stY:step:yLength;
   
   stX = mod(xLength-band , step);
   if(stX == 0)
       stX = step;
   end
   xLen = stX:step:xLength-band;
   
   stY = yLength-band+mod(band , step);
   yHeiA = stY:step:yLength;
   
   stX = xLength-band+mod(band , step);
   xHeiA = stX:step:xLength;
   
   %yLen = 1:step:yLength;
   %xLen = 1:step:floor(3*yLength/4);
   %yHeiA = floor(3*yLength/4):step:yLength;
   %xHeiA = floor(3*xLength/4):step:xLength;
end