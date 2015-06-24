function [uBe_i,PBe_i] = GetApproximationForInnerBoundary(x,y,c,boundaryU)

   
   xLength=length(x(ij:end));  yLength=length(y(lo:end));
 [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength);
%    yLen = 1:step:(floor((yLength-3*step - 1)/step)*step);
%    xLen = 1:step:(floor(xLength/step)*step);
%    yHeiA = yLength-3*step:step:yLength;
%    xHeiA = xLength-3*step:step:xLength;

   uBe_i=[boundaryU(xHeiA,yLen)  boundaryU(xLen,yHeiA)'];
   %cge_i=[cg(end-insb:end-8,1:msy)   cg(1:msx,end-insb:end-8)'];
   PBe_i = bt*(1-c^2)*uBe_i;
   
end