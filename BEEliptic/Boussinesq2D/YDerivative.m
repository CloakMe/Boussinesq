function zeroMatrix=YDerivative(M,zeroMatrix,augLeftFunctionPoints, augRightFunctionPoints,finiteDiff)   %dh operator if gm = 0
%{
%coeff table for 8th order finite diff
       %        -4         -3      -2     -1         0         1        2       3        4  
       %      -1/560     8/315   -1/5    8/5     -205/72      8/5     -1/5    8/315   -1/560
%}
    len = length(finiteDiff);
    mid = (len+1)/2;
    
%     j=1;
%     vdah(:,j) = uyst(:,1:4)*fofd(1:4)' + M(:,j:j+4)*fofd(5:9)';
%     j=2;
%     vdah(:,j) = uyst(:,2:4)*fofd(1:3)' + M(:,j-1:j+4)*fofd(4:9)';
%     j=3;
%     vdah(:,j) = uyst(:,3:4)*fofd(1:2)' + M(:,j-2:j+4)*fofd(3:9)';
%     j=4;
%     vdah(:,j) =  uyst(:,4)*fofd(1) + M(:,j-3:j+4)*fofd(2:9)';
    
    for j = 1:mid-1
        zeroMatrix(:,j) = augLeftFunctionPoints(:,j:mid-1)*finiteDiff(1:mid-j)' +...
            M(:,1:mid-1+j)*finiteDiff(mid+1-j:end)';
    end
    for j=mid:size(M,2)-mid+1
        zeroMatrix(:,j) = M(:,j-mid+1:j+mid-1)*finiteDiff';
    end

    for j = 0:mid-2
        zeroMatrix(:,end-j) =  M(:,end-mid+1-j:end)*finiteDiff(1:end-mid+1+j)' +...
            augRightFunctionPoints(:,1:mid-1-j)*finiteDiff(end-mid+2+j:end)';
    end

end