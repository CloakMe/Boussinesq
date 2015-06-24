function zeroMatrix=XDerivative(M,zeroMatrix,augFunctionPointsLeft,augFunctionPointsRigth,finiteDiff)
    
    zeroMatrix = YDerivative(M',zeroMatrix',augFunctionPointsLeft',augFunctionPointsRigth',finiteDiff)';
end