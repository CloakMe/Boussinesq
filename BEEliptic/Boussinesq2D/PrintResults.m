function PrintResults(solutionNorms,mu)
    
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
    
    fn = fieldnames(mu);
	for k=1:numel(fn)
	    if( isnumeric(mystruct.(fn{k})) )
	        fprintf('%s = %.4e \n', fn{1}, mu.(fn{1}));
	    end
	end

end