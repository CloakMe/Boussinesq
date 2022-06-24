function zeroMatrix=YDerivativeEvenFunctions(M,zeroMatrix,augFunctionPoints,finiteDiff)   %dh operator if gm = 0
%{
%coeff table for 8th order finite diff
       %        -4         -3      -2     -1         0         1        2       3        4  
       %      -1/560     8/315   -1/5    8/5     -205/72      8/5     -1/5    8/315   -1/560
%}
    len = length(finiteDiff);
    mid = (len+1)/2;

    for j = 1:mid-1
        zeroMatrix(:,j) =  M(:,1:mid-1+j)*finiteDiff(mid+1-j:end)'  + M(:,2:mid+1-j)*finiteDiff(mid+j:end)';
    end
    for j=mid:size(M,2)-mid+1
        zeroMatrix(:,j) = M(:,j-mid+1:j+mid-1)*finiteDiff';
    end
    
    h_square = sum(GetFiniteDifferenceCoeff(-(mid-1):(mid-1),2)'./finiteDiff)/len;
    if(sum(sum(augFunctionPoints)) == 0 && len == 5)
        c_i = GetFiniteDifferenceCoeff(-len:0,2)';
        d_i = GetFiniteDifferenceCoeff(-len+1:1,2)';
        d0_c0 = d_i(1)/c_i(1);
        customFiniteDiff = (d_i(2:end-1)-d0_c0*c_i(2:end-1))/h_square;
        j = 0;
        zeroMatrix(:,end-j) =  M(:,end-len+2:end)*customFiniteDiff';
        
        j = 1;
        zeroMatrix(:,end-j) =  M(:,end-mid+1-j:end)*finiteDiff(1:end-mid+1+j)';
        
    elseif(sum(sum(augFunctionPoints)) == 0 && len == 7)
        c_i = GetFiniteDifferenceCoeff(-len:0,2)';
        d_i = GetFiniteDifferenceCoeff(-len+1:1,2)';
        d0_c0 = d_i(1)/c_i(1);
        customFiniteDiff = (d_i(2:end-1)-d0_c0*c_i(2:end-1))/h_square;
        j = 0;
        zeroMatrix(:,end-j) =  M(:,end-len+2:end)*customFiniteDiff';
        
        c_i = GetFiniteDifferenceCoeff(-len:0,2)';
        d_i = GetFiniteDifferenceCoeff(-len+2:2,2)';
        d0_c0 = d_i(1)/c_i(1);
        customFiniteDiff = (d_i(2:end-1)-d0_c0*c_i(2:end-1))/h_square;
        j = 1;
        zeroMatrix(:,end-j) =  M(:,end-len+2:end)*customFiniteDiff';
                
        j = 2;
        zeroMatrix(:,end-j) =  M(:,end-mid+1-j:end)*finiteDiff(1:end-mid+1+j)';
        
    else
        for j = 0:mid-2
            zeroMatrix(:,end-j) =  M(:,end-mid+1-j:end)*finiteDiff(1:end-mid+1+j)' +  augFunctionPoints(:,1:mid-1-j)*finiteDiff(end-mid+2+j:end)';
        end
    end
end

