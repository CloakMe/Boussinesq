function vdah=DeltaEvenFunctions2D(M,vdah,socfd)  
    vdah= YDerivativeEvenFunctions2D(M,vdah,socfd) + XDerivativeEvenFunctions2D(M,vdah,socfd);
end