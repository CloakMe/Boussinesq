function [boundaryHit] = IsBoundaryHit(checkBnd,crrntResidual,residualInfNorm,subCounter,finiteDifference)
%    minResidual = min(minResidual,residualInfNorm(subCounter));
   if(checkBnd == 0)
       boundaryHit = 0;
       return;
   end
   mid = floor(length(finiteDifference)/2);
   z_B = max(max(abs(crrntResidual(end-mid:end,1:end))));
   z_B = max(z_B, max(max(abs(crrntResidual(1:end,end-mid:end)))));
   if(residualInfNorm(subCounter)<z_B/2 && z_B < 0.001)
%        figure(3)
%        mesh(x(zeroX:end),y(zeroY:end),crrntResidual');
       warning('Boundary Condition Cant Improve! Stopping!');
       boundaryHit = 1;
       return;
   end
   boundaryHit = 0;
end