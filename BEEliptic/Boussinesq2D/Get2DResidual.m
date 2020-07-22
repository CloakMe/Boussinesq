function crrntResidual = Get2DResidual(al,bt,c,th,iterCounter,U,Usquare,deltaU,...
    zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,secondDerivative)

    dyy_yk = YDerivativeEvenFunctions(U,zeroMatrix,outerTopBoundaryF,secondDerivative);
    crrntResidual =  DeltaEvenFunctions(bt*U + al*bt*th(iterCounter)*Usquare - deltaU +...
       bt*c^2*dyy_yk ,zeroMatrix,bt*outerRigthBoundaryF,bt*outerTopBoundaryF,secondDerivative) -...
       bt*c^2*dyy_yk;
   
end
           