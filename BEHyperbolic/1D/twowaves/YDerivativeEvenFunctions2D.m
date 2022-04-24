function zeroMatrix=YDerivativeEvenFunctions2D(M,zeroMatrix,finiteDiff)   %dh operator if gm = 0
%{
%coeff table for 8th order finite diff
       %        -4         -3      -2     -1         0         1        2       3        4  
       %      -1/560     8/315   -1/5    8/5     -205/72      8/5     -1/5    8/315   -1/560
%}
    len = length(finiteDiff);
    mid = (len+1)/2;

    for j = 1:mid-1
        zeroMatrix(:,j) =  M(:,1:mid-1+j)*finiteDiff(mid+1-j:end)'  + M(:,2:mid+1-j)*finiteDiff(mid-j:-1:1)';
    end
    for j=mid:size(M,2)-mid+1
        zeroMatrix(:,j) = M(:,j-mid+1:j+mid-1)*finiteDiff';
    end
    for j = 0:mid-2
        zeroMatrix(:,end-j) =  M(:,end-mid+1-j:end)*finiteDiff(1:mid+j)' +  M(:,end-mid+1+j:end-1)*finiteDiff(mid+1+j:end)';
    end
end

