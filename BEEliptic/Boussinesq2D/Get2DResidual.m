function crrntResidual = Get2DResidual(al,bt,c,th,iterCounter,U,Usquare,deltaU,...
    zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,secondDerivative, prmtrs, domUtils)
        
    if(nargin == 14 && prmtrs.useZeroBoundary == 2)
        dyy_yk = domUtils.YDerivativeEvenFunZeroBnd(U)/prmtrs.h^2;
        crrntResidual =  domUtils.DeltaHEvenFunZeroBnd(bt*U + al*bt*th(iterCounter)*Usquare - deltaU +...
            bt*c^2*dyy_yk )/prmtrs.h^2 -...
            bt*c^2*dyy_yk;
    else
        dyy_yk = YDerivativeEvenFunctions(U,zeroMatrix,outerTopBoundaryF,secondDerivative);
        crrntResidual =  DeltaEvenFunctions(bt*U + al*bt*th(iterCounter)*Usquare - deltaU +...
            bt*c^2*dyy_yk ,zeroMatrix,bt*outerRigthBoundaryF,bt*outerTopBoundaryF,secondDerivative) -...
            bt*c^2*dyy_yk;
    end
   
end
           