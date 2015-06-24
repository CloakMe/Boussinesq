function zeroMatrix=YDerivativeEvenFunctions2(M,zeroMatrix,augFunctionPoints,finiteDiff)   %dh operator if gm = 0
%{
%coeff table for 8th order finite diff
       %        -4         -3      -2     -1         0         1        2       3        4  
       %      -1/560     8/315   -1/5    8/5     -205/72      8/5     -1/5    8/315   -1/560
%}
    len = length(finiteDiff);
    mid = (len+1)/2;

    for j = 0:mid-2
        zeroMatrix(:,j+1) =  M(:,1:mid+j)*finiteDiff(mid-j:end)'  + M(:,2:mid-j)*finiteDiff(mid+j+1:end)';
        zeroMatrix(:,end-j) =  M(:,end-mid+1-j:end)*finiteDiff(1:end-mid+1+j)' +...
            augFunctionPoints(:,1:mid-1-j)*finiteDiff(end-mid+2+j:end)';
    end
    mat3 = zeros([size(M,1) (size(M,1)-len+1) len]);
    for j=0:len-1
        mat3(:,:,j+1) =finiteDiff(j+1).* M(:,1+j:end-2*mid+2+j);
    end
       zeroMatrix(:,mid:end-mid+1) = sum(mat3,3);
%     for j=mid:size(M,2)-mid+1
%         zeroMatrix(:,j) = M(:,j-mid+1:j+mid-1)*finiteDiff';
%     end

%     for j = 0:mid-2
%         
%     end

end
% size(zeroMatrix(:,mid:end-mid+1))
% size(sum(mat3,3))