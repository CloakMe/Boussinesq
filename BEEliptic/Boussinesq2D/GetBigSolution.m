function [bigU,bigUTimeDerivative] = GetBigSolution(x,y,c,bt,thetaEnd,U,bigZeroMatrix,firstDerivative,outerTopBoundaryF)

    bigOuterTopBoundaryF = fliplr([flipud(outerTopBoundaryF(2:end,:)); outerTopBoundaryF]);
    bigOuterBottomBoundaryF = [flipud(outerTopBoundaryF(2:end,:)); outerTopBoundaryF];
    
    [zeroX,zeroY] = GetZeroNodes(x,y);
    bigU = thetaEnd*transf2qD(U,x,y,zeroX,zeroY);
    
    bigUTimeDerivative = -c * YDerivative(bigU,bigZeroMatrix,thetaEnd*bigOuterTopBoundaryF,...
        thetaEnd*bigOuterBottomBoundaryF,firstDerivative); 
end