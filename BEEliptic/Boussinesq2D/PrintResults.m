function PrintResults(solutionNorms,c1,c2)
    
    fprintf('||Residual||_Inf at the start = %.4e \n', solutionNorms.residualInfNorm(1));
    fprintf('At the end:\n');
    fprintf('||Residual||_Inf = %.4e \n', solutionNorms.residualInfNorm(end));
    fprintf('||Residual||_L2 = %.4e \n\n', solutionNorms.residualL2Norm);

    
    fprintf('At the end:\n');
    fprintf('||P-Pup||_Inf = %.4e \n', solutionNorms.PvsPupInfNorm);
    fprintf('||U-Uup||_Inf = %.4e \n', solutionNorms.UvsUpInfNorm(end));

    fprintf('||P-Pup||_L2 = %.4e \n', solutionNorms.PvsPupL2Norm);
    fprintf('||U-Uup||_L2 = %.4e \n', solutionNorms.UvsUupL2Norm);

    fprintf('||PBndAprox - PBnd||_L2 = %.4e \n', solutionNorms.BoundaryFunctionPvsPL2Norm);
    fprintf('||UBndAprox - UBnd||_L2 = %.4e \n\n', solutionNorms.boundaryFunctionUvsUL2Norm);
    fprintf('c1 = %.4e \n', c1);
    fprintf('c2 = %.4e \n', c2);

end