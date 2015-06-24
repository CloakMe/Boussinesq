function [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step)
    
% The approximation boundary function (derived from the simplified BPE)
% and the current solution (with respect
% to the current iteration) are to be compared on the following node numbers:
% (xHeiA,yLen) and (xLen,yHeiA)
  
   yLen = 1:step:(floor((yLength-3*step - 1)/step)*step);
   xLen = 1:step:(floor(xLength/step)*step);
   yHeiA = yLength-3*step:step:yLength;
   xHeiA = xLength-3*step:step:xLength;

end