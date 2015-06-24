function zeroMatrix=XDer(M,finiteDiff)   %dh operator if gm = 0
    zeroMatrix=YDer(M',finiteDiff)';
end