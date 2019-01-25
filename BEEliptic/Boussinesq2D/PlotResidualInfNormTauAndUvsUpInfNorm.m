function PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl)
figure(1)
len = length(tauVector);

loglog(1:10:10*length(solutionNorms.residualInfNorm),solutionNorms.residualInfNorm,...
    (10:10:len),tauVector(10:10:len),'g',...
    (1:10:len),solutionNorms.UvsUpInfNorm(1:10:len),'r');%,...
    %(21:10:len+1),angl(3:end)-180,'c');
title('Residual, tau and UvsUupInfNorm');
legend('|Residual|_{Inf} Norm','tau','|U-Up|_{Inf} Norm')
end
%{
figure(5);
len2 = length(rsdlInf(1:subIter));
plot(1:len2,rsdlInf(1:subIter));
%}
%figure(14)
%plot(len-10000:len, tauVector(len-10000:len))