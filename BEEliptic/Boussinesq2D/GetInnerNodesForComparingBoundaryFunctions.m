function [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step)
    
% The approximation boundary function (derived from the simplified BPE)
% and the current solution (with respect
% to the current iteration) are to be compared on the following node numbers:
% (xHeiA,yLen) and (xLen,yHeiA)

   yLen = 1:step:yLength;
   xLen = 1:step:floor(3*yLength/4);
   yHeiA = floor(3*yLength/4):step:yLength;
   xHeiA = floor(3*xLength/4):step:xLength;
end