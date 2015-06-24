function [rightBoundaryU,topBoundaryU] = GetApproximationForOuterBoundary(boundaryU)

   rightBoundaryU=boundaryU(1:end-4,end-3:end); topBoundaryU=boundaryU(end-3:end,1:end-4); 
  
end