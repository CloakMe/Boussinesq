%return;
clear;clc;
addpath('..\Boussinesq2D');
tic
x_st = -20.0;
x_end = 20.0;
y_st = -20.0;
y_end = 20.0;

h = 0.4;
x=x_st:h:x_end; 
y=y_st:h:y_end;

sx = (length(x)+1)/2;
fprintf('x size = %d\n', sx);
sy = (length(y)+1)/2;
fprintf('y size = %d\n', sy);

al = 1;%99979 izb
bt1 = 2;bt2 = 1; bt = bt1/bt2;
c = 1.25; 
iterMax = 9000000;
%eps = 1/max(y_end^6,((1-c^2)*x_end^2)^3);
eps = 5.0e-04;%5.0e-09;
tau = getTau(h,x_end,y_end)/200;
plotResidual  = 0;
plotBoundary  = 0;
plotAssympt   = 0;
prmtrs = struct('h',{h},'tau',{tau},'iterMax',{iterMax},'eps',{eps},...
    'plotResidual',{plotResidual},'plotBoundary',{plotBoundary},'plotAssympt',{plotAssympt});

firstDerivative = GetFiniteDifferenceCoeff([-1,0,1],1)'/h;
secondDerivative = GetFiniteDifferenceCoeff([-1,0,1],2)'/h^2;
derivative = struct('first',{firstDerivative},'second',{secondDerivative});

[X,Y]=Domain(x,y);
%IC = 3/2 * (c^2 - 1) * sech(sqrt(bt1 * (c ^ 2 - 1) / (bt1 * c ^ 2 - bt2)) * X / 2) .^ 2 / al +...
%     3/2 * (c^2 - 1) * sech(sqrt(bt1 * (c ^ 2 - 1) / (bt1 * c ^ 2 - bt2)) * Y / 2) .^ 2 / al;
k = 1.25; % k1 = -.3; k2 = .3; 
%b1 := 4; b2 := 4;
a1 = .5; a2 = .5; a12 = .25;
IC = GetApproximateSolution(X,Y,k,a1,a2,a12);
figure(10)
mesh(x,y,IC');
xlabel('s');    ylabel('r');
title('start solution');
return;
[zeroX,zeroY]=GetZeroNodes(x,y);
th = abs(IC(zeroX,zeroY));
IC = IC/th;
%U = -5*tanh(x/2).^2+5;

[Phi,Psi,thetaVector,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_1d_v2(IC,x,y,prmtrs,bt1,bt2,al,c,th,derivative);

        figure(20)
        plot(x,Phi(1,:));
        
fprintf('elapsed time = %d \n', toc);

figure(1)
mesh(x, y, Phi);
title('end solution');
figure(2)
mesh(x, y, Psi);
title('end supplementary function');
figure(3);
resEnd = ceil(length(tauVector)/10);
plot(1:resEnd, solutionNorms.residualInfNorm(1:resEnd), 'b', 2:length(tauVector), solutionNorms.UvsUpInfNorm(2:end), 'r' );
title('residualInfNorm');
figure(4);
plot(1:length(thetaVector), thetaVector, 'g' );
title('theta');
%   if(sw_div == 1)
%         return;
%   end
%   save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
%       '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
  


return;

% Continue from lasth iteration:
lastTheta=theta(end); lastU=U; lastP = P;  last_tau = tauVector(end); 

[PhiUp,PsiUp,thetaVector,solutionNorms,tauVector,angl,sw_div] =...
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
