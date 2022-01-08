%return;
clear;clc;
addpath('..\Boussinesq2D');
tic
x_st = 0.0;
x_end = 45.0;
y_st = 0.0;
y_end = 45.0;

h = 0.1;
x=x_st:h:x_end; 
y=y_st:h:y_end;

sx = (length(x)+1)/2;
fprintf('x size = %d\n', sx);
sy = (length(y)+1)/2;
fprintf('y size = %d\n', sy);

al = 1;%99979 izb
bt1 = 3;bt2 = 1; bt = bt1/bt2;
c = 2; 
iterMax = 9000000;
%eps = 1/max(y_end^6,((1-c^2)*x_end^2)^3);
eps = 5.0e-04;%5.0e-09;
tau = getTau(h,x_end,y_end)/100;

prmtrs = struct('h',{h},'tau',{tau},'iterMax',{iterMax},'eps',{eps});

firstDerivative = GetFiniteDifferenceCoeff([-2,-1,0,1,2],1)'/h;
secondDerivative = GetFiniteDifferenceCoeff([-2,-1,0,1,2],2)'/h^2;
derivative = struct('first',{firstDerivative},'second',{secondDerivative});

[X,Y]=Domain(x,y);
IC = 3/2 * (c^2 - 1) * sech(sqrt(bt1 * (c ^ 2 - 1) / (bt1 * c ^ 2 - bt2)) * X / 2) .^ 2 / al +...
     3/2 * (c^2 - 1) * sech(sqrt(bt1 * (c ^ 2 - 1) / (bt1 * c ^ 2 - bt2)) * Y / 2) .^ 2 / al;
figure(10)
mesh(x,y,IC');
xlabel('s');    ylabel('r');
title('start solution');


th = abs(IC(1,1));
IC = IC/th;
%U = -5*tanh(x/2).^2+5;

[Fup,Uup,residualInfNorm,thetaVector,tauVector]=...
    sol_ch_1d_v2(IC,x,y,prmtrs,bt1,bt2,al,c,th,derivative);

figure(1)
plot(x,Uup);
figure(2)
plot(x,Fup);
figure(3);
plot(1:length(residualInfNorm),residualInfNorm);

%   if(sw_div == 1)
%         return;
%   end
%   save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
%       '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
  
fprintf('elapsed time = %d \n', toc);

return;

% Continue from lasth iteration:
lastTheta=theta(end); lastU=U; lastP = P;  last_tau = tauVector(end); 

[bigU,bigUTimeDerivative,P,U,thetaCont,mu,solutionNormsCont,tauVecCont,anglCont,sw_div] =...
       sol_ch_v8(lastU,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastP);
if(sw_div == 1)
    return;
end
tauVector = [tauVector tauVecCont];
theta = [ theta thetaCont];
angl = [angl anglCont];

solutionNorms.residualInfNorm = [solutionNorms.residualInfNorm solutionNormsCont.residualInfNorm];
solutionNorms.UvsUpInfNorm = [solutionNorms.UvsUpInfNorm solutionNormsCont.UvsUpInfNorm];
solutionNorms.residualL2Norm = solutionNormsCont.residualL2Norm;
solutionNorms.UvsUupL2Norm = solutionNormsCont.UvsUupL2Norm;
solutionNorms.PvsPupInfNorm = solutionNormsCont.PvsPupInfNorm;
solutionNorms.PvsPupL2Norm = solutionNormsCont.PvsPupL2Norm;
solutionNorms.boundaryFunctionUvsUL2Norm = solutionNormsCont.boundaryFunctionUvsUL2Norm;
solutionNorms.BoundaryFunctionPvsPL2Norm = solutionNormsCont.BoundaryFunctionPvsPL2Norm;

save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
      '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
PrintResults(solutionNorms,mu);
PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c);
return;

figure(1)
mesh(x(zeroX:end),y(zeroY:end),UTilda');
figure(2)
mesh(x(zeroX:end),y(zeroY:end),Uup');
figure(3)
mesh(x(zeroX:end),y(zeroY:end),Uup1');        
