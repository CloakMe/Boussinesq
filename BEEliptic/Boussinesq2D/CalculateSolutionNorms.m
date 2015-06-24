function [solutionNorms] = CalculateSolutionNorms(U,Uup,P,Pup,UvsUupInfNorm,...
    crrntResidual,residualInfNorm,innerBoundaryF,innerBoundaryPF,subCounter,thetaEnd,step,h,c1,c2)
    
    %Already multiplied by thetaEnd
    residualInfNorm = residualInfNorm(1:subCounter);
    
    [boundaryFunctionUvsUL2Norm,BoundaryFunctionPvsPL2Norm]=...
    CompareBoundaryFunctionAndIterSolution(Uup,Pup,c1*innerBoundaryF,c2*innerBoundaryPF,step,h);
    
    residualL2Norm = thetaEnd*h*norm(crrntResidual,2);
    PvsPupInfNorm = thetaEnd*max(max(abs(Pup-P)));
    PvsPupL2Norm = thetaEnd*h*norm(Pup(:)-P(:),2);
    UvsUupL2Norm = thetaEnd*h*norm(Uup(:)-U(:),2);
    
    solutionNorms = struct('residualInfNorm',{residualInfNorm},'residualL2Norm',{residualL2Norm},...
     'UvsUpInfNorm',{UvsUupInfNorm},'UvsUupL2Norm',{UvsUupL2Norm},'PvsPupInfNorm',{PvsPupInfNorm},...
     'PvsPupL2Norm',{PvsPupL2Norm},'boundaryFunctionUvsUL2Norm',{thetaEnd*boundaryFunctionUvsUL2Norm},...
     'BoundaryFunctionPvsPL2Norm',{thetaEnd*BoundaryFunctionPvsPL2Norm});
end