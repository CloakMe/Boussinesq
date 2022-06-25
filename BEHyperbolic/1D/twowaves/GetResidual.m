function crrntResidual = GetResidual(type, bt, btExt, c, Phi, coefPhiSquare, deltaPhi, zeroMatrix, derivative)
    if(strcmp(type,'pq') == 1)
        Psi = (-1) * deltaPhi -...
            2*(+1) * XDerivative( YDerivative( Phi, zeroMatrix, derivative.firstT), zeroMatrix, derivative.firstX) + coefPhiSquare;  

        crrntResidual = bt*(c^2 - 1) * deltaPhi - XDerivative(Psi, zeroMatrix,derivative.second) - YDerivative(Psi, zeroMatrix,derivative.secondT) -...
            XDerivative( YDerivative( 2*bt*(c^2 + 1) * Phi + 2 * Psi, zeroMatrix, derivative.firstT), zeroMatrix, derivative.firstX);
    elseif(strcmp(type, 'xt') == 1)

        crrntResidual = - bt * XDerivative(Phi, zeroMatrix,derivative.secondX) + bt * YDerivative(Phi, zeroMatrix,derivative.secondT) -...
            XDerivative( coefPhiSquare*bt*Phi.^2, zeroMatrix,derivative.secondX )+ ...
            XDerivative(XDerivative(Phi, zeroMatrix,derivative.secondX), zeroMatrix,derivative.secondX);
    elseif(strcmp(type, 'xte') == 1)

        crrntResidual = - bt * XDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondX) + bt * YDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondT) -...
            XDerivativeEvenFunctions2D( coefPhiSquare*bt*Phi.^2, zeroMatrix,derivative.secondX ) + ...
            XDerivativeEvenFunctions2D(XDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondX), zeroMatrix,derivative.secondX) - ...
            btExt * bt * XDerivativeEvenFunctions2D(YDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondT), zeroMatrix,derivative.secondX);
   elseif(strcmp(type, 'org') == 1)

        crrntResidual = - bt * XDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondX) + bt * YDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondT) -...
            XDerivativeEvenFunctions2D( coefPhiSquare*bt*Phi.^2, zeroMatrix,derivative.secondX ) + ...
            XDerivativeEvenFunctions2D(XDerivativeEvenFunctions2D(Phi, zeroMatrix,derivative.secondX), zeroMatrix,derivative.secondX);
    else
        error('unknown type');
    end
end