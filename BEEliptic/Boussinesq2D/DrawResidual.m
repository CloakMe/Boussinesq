function crrntResidual = DrawResidual(x,y,al,bt,c,th,c1,U,secondDerivative)

    zeroMatrix = zeros(size(U));
    Usquare = U .^2;
    deltaU = DeltaEvenFunctions(U, zeroMatrix, zeroMatrix, zeroMatrix, secondDerivative);
    Uyy = YDerivativeEvenFunctions(U,zeroMatrix,zeroMatrix,secondDerivative);
    crrntResidual =  DeltaEvenFunctions(bt*U + al*bt*th*Usquare - deltaU +...
       bt*c^2*Uyy ,zeroMatrix,zeroMatrix,zeroMatrix,secondDerivative) -...
       bt*c^2*Uyy;
   
   figure(3)
    mesh(x(1:end-3),y(1:end-3),crrntResidual((1:end-3),(1:end-3))');
    xlabel('x');    ylabel('y');
    view(-13,30);
    title('Residual');
   
end