function crrntResidual = GetResidual(bt, c, Phi, coefPhiSquare, deltaPhi, zeroMatrix, derivative)

	Psi = bt*(c^2-1) * DeltaEvenFunctions2D(Phi, zeroMatrix,derivative.second) -...
             2*bt*(c^2+1) * XDerivativeEvenFunctions2D( YDerivativeEvenFunctions2D( Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first);  

	crrntResidual = bt*(c^2 - 1) * deltaPhi - DeltaEvenFunctions2D(Psi, zeroMatrix,derivative.second) + coefPhiSquare -...
        XDerivativeEvenFunctions2D( YDerivativeEvenFunctions2D( 2*bt*(c^2 + 1) * Phi + 2 * Psi, zeroMatrix, derivative.first), zeroMatrix, derivative.first);
end