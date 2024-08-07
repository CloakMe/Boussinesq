%return;
clear;clc;
%addpath('..\2D\OOPDesign');
tic
x_st =  -40.0;
x_end = 40.0;
y_st =  0.0;
y_end = 0.0;

h = 0.2;
tTau = 0.2;
x=x_st:h:x_end; 
y=y_st:tTau:y_end;

sx = (length(x)+1)/2;
fprintf('x size = %d\n', sx);
sy = (length(y)+1)/2;
fprintf('y size = %d\n', sy);

al = 3;%99979 izb
bt1 = 1;bt2 = 1; bt = bt1/bt2;
btExt = 0.1;
c = 0.2862; 
iterMax = 9000000;
%eps = 1/max(y_end^6,((1-c^2)*x_end^2)^3);
eps = 5.0e-04;%5.0e-09;
tau = 1; %getTau(h,x_end,y_end)/500;
plotResidual  = 1;
plotBoundary  = 0;
plotAssympt   = 0;
type = 'pq';
prmtrs = struct('h',{h},'tTau',{tTau},'tau',{tau},'iterMax',{iterMax},'eps',{eps}, 'type', {type}, ...
    'plotResidual',{plotResidual},'plotBoundary',{plotBoundary},'plotAssympt',{plotAssympt},'btExt',{btExt});

firstXDerivative = GetFiniteDifferenceCoeff([-1,0,1],1)'/h;
secondXDerivative = GetFiniteDifferenceCoeff([-1,0,1],2)'/h^2;
firstTDerivative = GetFiniteDifferenceCoeff([-1,0,1],1)'/tTau;
secondTDerivative = GetFiniteDifferenceCoeff([-1,0,1],2)'/tTau^2;
derivative = struct('firstX',{firstXDerivative},'secondX',{secondXDerivative},'firstT',{firstTDerivative},'secondT',{secondTDerivative});
% firstDerivative = firstDerivative
% secondDerivative = secondDerivative
% return;
[X,Y]=Domain(x,y);
%IC = 3/20 * (c^2 - 1) * sech(sqrt(bt1 * (c ^ 2 - 1) / (bt1 * c ^ 2 - bt2)) * X / 2) .^ 2 / al +...
%    3/20 * (c^2 - 1) * sech(sqrt(bt1 * (c ^ 2 - 1) / (bt1 * c ^ 2 - bt2)) * Y / 2) .^ 2 / al;
k = sqrt(.5 - sqrt(1-4*c^2)/2); % k1 = -.3; k2 = .3; 
%b1 := 4; b2 := 4;
a1 = 0.5; a2 = 0.5; a12 = 0.25;
IC = 1*GetApproximateSolution('pq',X,Y,k,a1,a2,a12);
figure(1)
mesh(x,y,IC');
xlabel('x');    ylabel('t');
title('start solution');
return;
zeroMatrix = zeros(size(IC));

%[zeroX,zeroY]=GetZeroNodes(x,y);
% th = abs(IC(zeroX,zeroY));
% IC = IC/th;
th = 1;
crrntResidual = th*GetResidual(type, bt, btExt, c, IC, al*th, 0, zeroMatrix, derivative);
figure(1)
%mesh(x,y,crrntResidual');
mesh(x(3:end-2), y(5:end-4), crrntResidual(3:end-2,5:end-4)');
xlabel('x');    ylabel('t');
title('start residual');
% return;

[Phi,Psi,thetaVector,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_1d_v6(IC,x,y,prmtrs,bt1,bt2,al,c,th,derivative);

fprintf('elapsed time = %d \n', toc);

 figure(10)
 %mesh(x, y, Phi');
 mesh(x(1:end), y(2:end-1), Phi(1:end,2:end-1)');
 xlabel('x');    ylabel('t');
 title('end solution');
% figure(11)
% mesh(x, y, (Phi-IC)');
% title('diff end solution - IC');
% figure(12)
% %mesh(x, y, Psi');
% mesh(x(1:end), y(3:end-2), Psi(1:end,3:end-2)')
% xlabel('x');    ylabel('t');
% title('end supplementary function');
figure(2)
%deltaPhi = XDerivative(Phi, zeroMatrix,derivative.second) + YDerivative(Phi, zeroMatrix,derivative.second);
crrntResidual = th*GetResidual(type, bt, btExt, c, Phi, al*th, 0, zeroMatrix, derivative);
mesh(x(3:end-2), y(5:end-4), crrntResidual(3:end-2,5:end-4)');
title('end resiudal');
return;
figure(5);
resEnd = ceil(length(tauVector)/10);
plot(1:resEnd, solutionNorms.residualInfNorm(1:resEnd), 'b', 1:resEnd, solutionNorms.UvsUpInfNorm(1:10:end), 'r' );
title('residualInfNorm');
figure(6);
plot(1:length(thetaVector), thetaVector, 'g' );
title('theta');
%   if(sw_div == 1)
%         return;
%   end
%   save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
%       '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
  


return;

% Continue from lasth iteration:
lastTheta=thetaVector(end); lastPhi=Phi; lastPsi = Psi;  last_tau = tauVector(end); 

[PhiUp,PsiUp,thetaVector,solutionNorms,tauVector,angl,sw_div] =...
	sol_ch_1d_v2(lastPhi,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastPsi);
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

figure(1)
mesh(x(1:end), y(15:end-15), Phi(1:end,15:end-15)');
xlabel('x');    ylabel('t');
title('end solution');
figure(2)
mesh(x(15:end-15), y(15:end-15), Psi(15:end-15,15:end-15));
title('end supplementary function');
figure(3);
resEnd = ceil(length(tauVector)/10);
plot(1:resEnd, solutionNorms.residualInfNorm(1:resEnd), 'b', 2:length(tauVector), solutionNorms.UvsUpInfNorm(2:end), 'r' );
title('residualInfNorm');
figure(5);
plot(1:length(thetaVector), thetaVector, 'g' );
title('theta');

%save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
%      '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
%PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
%PrintResults(solutionNorms,mu);
%PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c);
return;

figure(1)
mesh(x(zeroX:end),y(zeroY:end),UTilda');
figure(2)
mesh(x(zeroX:end),y(zeroY:end),Uup');
figure(3)
mesh(x(zeroX:end),y(zeroY:end),Uup1');        
