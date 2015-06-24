function zeroMatrix=DeltaEvenFunctions(M,zeroMatrix,augFunctionPointsX,augFunctionPointsY,finiteDiff)
        
    zeroMatrix = YDerivativeEvenFunctions(M',zeroMatrix',augFunctionPointsX',finiteDiff)' +...
        YDerivativeEvenFunctions(M,zeroMatrix,augFunctionPointsY,finiteDiff);
end