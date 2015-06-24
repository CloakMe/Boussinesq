function BndStats(U,P,solutionNorms,c1,c2)

    fprintf('U(xEnd,0) = %.4e \n', U(end,1));
    fprintf('U(0,yEnd) = %.4e \n', U(1,end));
    
    fprintf('P(xEnd,0) = %.4e \n', P(end,1));
    fprintf('P(0,yEnd) = %.4e \n\n', P(1,end));
    
    fprintf('||UBoundary - U_L2||_L2 = %.4e \n', solutionNorms.boundaryFunctionUvsUL2Norm);
    fprintf('||PBoundary - P_L2||_L2 = %.4e \n', solutionNorms.BoundaryFunctionPvsPL2Norm);
        
    fprintf('mu1 = %.4e \n', c1);
    fprintf('mu2 = %.4e \n', c2);

end