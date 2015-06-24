function crrntResidual = Get2DResidual(al,bt,c,th,c1,iterCounter,U,Usquare,deltaU,...
    zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,secondDerivative)

    dyy_yk = YDerivativeEvenFunctions(U,zeroMatrix,c1*outerTopBoundaryF,secondDerivative);
    crrntResidual =  DeltaEvenFunctions(bt*U + al*bt*th(iterCounter)*Usquare - deltaU +...
       bt*c^2*dyy_yk ,zeroMatrix,bt*c1*outerRigthBoundaryF,bt*c1*outerTopBoundaryF,secondDerivative) -...
       bt*c^2*dyy_yk;
   
end
           