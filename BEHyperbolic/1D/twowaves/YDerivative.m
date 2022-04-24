function zeroMatrix=YDerivative(M,zeroMatrix,finiteDifferences)   %dh operator if gm = 0
%{
%coeff table for 8th order finite diff
       %        -4         -3      -2     -1         0         1        2       3        4  
       %      -1/560     8/315   -1/5    8/5     -205/72      8/5     -1/5    8/315   -1/560
%}
    fdSize = size(finiteDifferences);
    centralFiniteDiffPosition = (fdSize(1)+1)/2;
    midPoint = centralFiniteDiffPosition; %fdSize(2)/2;
    
    for j = 1:midPoint-1
        zeroMatrix(:,j) =  M(:,1:fdSize(2))*finiteDifferences(j,:)';
    end
    for j = midPoint:size(M,2)-midPoint+1
        zeroMatrix(:,j) = M(:,j-midPoint+1:j+midPoint-1)*finiteDifferences(centralFiniteDiffPosition,1:fdSize(1))';
    end
    for j = 0:midPoint-2
        zeroMatrix(:,end-j) =  M(:,end-fdSize(2)+1:end)*finiteDifferences(end-j,:)';
    end
end