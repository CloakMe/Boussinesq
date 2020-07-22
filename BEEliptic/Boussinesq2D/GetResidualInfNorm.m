function residualInfNorm=GetResidualInfNorm(al,bt,c,th,iterCounter,U,Usquare,deltaU,zeroMatrix,outerRigthBoundaryF,...
    outerTopBoundaryF,secondDerivative)

    mid = floor(length(secondDerivative)/2);
    crrntResidual = Get2DResidual(al,bt,c,th,iterCounter,U,Usquare,...
        deltaU,zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,secondDerivative);

    residualInfNorm=max(max(abs(crrntResidual(1:end-mid ,1:end-mid ))));
end