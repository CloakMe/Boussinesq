function zeroMatrix=DeltaEvenFunctions2(M,zeroMatrix,augFunctionPointsX,augFunctionPointsY,finiteDiff)
        
    zeroMatrix = YDerivativeEvenFunctions2(M',zeroMatrix',augFunctionPointsX',finiteDiff)' +...
        YDerivativeEvenFunctions2(M,zeroMatrix,augFunctionPointsY,finiteDiff);
end