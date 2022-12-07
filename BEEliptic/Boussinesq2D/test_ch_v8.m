%return;
addpath('..\..\Common');
clear;clc;
tic
x_st = -12.0;    y_st = -12.0;
x_end = 12.0;    y_end = 12.0;
x_st2 = -20.0;   y_st2 = -20.0;
x_end2 = 20.0;   y_end2 = 20.0;

compBox = struct('x_st',{x_st},'x_end',{x_end},'y_st',{y_st},'y_end',...
    {y_end},'x_st2',{x_st2},'x_end2',{x_end2},'y_st2',{y_st2},'y_end2',{y_end2});

UseExtendedDomain=1;

h = 0.5;
x=x_st2:h:x_end2; 
y=y_st2:h:y_end2; 
%tau = 0.00114425*8;% getTau(h,x_end,y_end)/20;
tau = getTau(h,x_end,y_end)/20;

sx = (length(x)+1)/2
sy = (length(y)+1)/2
   
   al = -1;%99979 izb
   bt1 = 3;bt2 = 1; bt = bt1/bt2;
   c = 0.45;
   iterMax = 9000000;
   %eps = 1/max(y_end^6,((1-c^2)*x_end^2)^3);
   eps = 5.0e-06;%5.0e-09;
   % IC_switch = 0 ->christov sech formula
   % IC_switch = 1 ->nat42 formula
   ICSwitch=0;   
   %useZeroBoundary:
   % 0 -> use points (one, two, three for p=2,4,6 ) outside the domain with nonzero boundary function
   % 1 -> use points (one, two, three for p=2,4,6 ) outside the domain with zero boundary function
   % 2 -> use non-symetric finite differences without points outside the domain   
   useZeroBoundary  = 0;
   plotResidual  = 0;
   plotBoundary  = 0;
   checkBoundary = 0;
   plotAssympt   = 0;
   % if '1' plots the Residual/Boundry
   % if '0' does not plot anything
   prmtrs = struct('h',{h},'tau',{tau},'iterMax',{iterMax},'eps',{eps},'ICSwitch',{ICSwitch},...
       'plotResidual',{plotResidual},'plotBoundary',{plotBoundary},'plotAssympt',{plotAssympt},...
       'checkBoundary',{checkBoundary}, 'useZeroBoundary', {useZeroBoundary});
   
   firstDerivative = GetFiniteDifferenceCoeff([-2,-1,0,1,2],1)'/h;
   secondDerivative = GetFiniteDifferenceCoeff([-2,-1,0,1,2],2)'/h^2;
   derivative = struct('first',{firstDerivative},'second',{secondDerivative});
   
  [bigU,bigUTimeDerivative,P,U,bigIC,solutionNorms,theta,mu,tauVector,angl,sw_div]=...
  PrepareICForInnerDomain(compBox,prmtrs,al,bt1,bt2,c,derivative);

  if(length(tauVector)<iterMax && UseExtendedDomain == 1 && size(bigUTimeDerivative,1)~=1)
     fprintf('\nLarge Domamin Calculations:\n\n');
     prmtrs.checkBoundary = 0;
     prmtrs.eps = 1.0e-13;
     prmtrs.plotResidual = 0;
     prmtrs.tau = tauVector(end);
     [bigU,bigUTimeDerivative,P,U,newBigIC,solutionNorms,theta,mu,tauVector,angl,sw_div] =...
     PrepareICForEnlargedDomain(bigU,compBox,prmtrs,al,bt1,bt2,c,mu,theta(end),derivative);
     x=x_st2:h:x_end2; y=y_st2:h:y_end2;
  end  
  if(sw_div == 1)
        return;
  end
  save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
      '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
  
fprintf('elapsed time = %d \n', toc);
PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
PrintResults(solutionNorms,mu);
PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c, 2);
return;

% Continue from lasth iteration:
lastTheta=theta(end); lastU=U; lastP = P;  last_tau = tauVector(end); 

if(prmtrs.useZeroBoundary == 2)    
    [bigU,bigUTimeDerivative,P,U,thetaCont,mu,solutionNormsCont,tauVecCont,anglCont,sw_div] =...
        sol_ch_v8ZeroBnd(lastU,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastP);
else
    [bigU,bigUTimeDerivative,P,U,thetaCont,mu,solutionNormsCont,tauVecCont,anglCont,sw_div] =...
        sol_ch_v8(lastU,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastP);
end
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
solutionNorms.UL2Norm = solutionNormsCont.UL2Norm;
solutionNorms.PL2Norm = solutionNormsCont.PL2Norm;
solutionNorms.boundaryFunctionUvsUInfNorm = solutionNormsCont.boundaryFunctionUvsUInfNorm;
solutionNorms.BoundaryFunctionPvsPInfNorm = solutionNormsCont.BoundaryFunctionPvsPInfNorm;
solutionNorms.UInfNorm = solutionNormsCont.UInfNorm;
solutionNorms.PInfNorm = solutionNormsCont.PInfNorm;

save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
      '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);
PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
PrintResults(solutionNorms,mu);
PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c);
return;

% Create smoother solution from existing save:
nh = h/2;
[bigU,bigUTimeDerivative,P,U,theta,c1,c2,solutionNorms,tauVector,angl,x,y,h] =...
    PreSolverFromExistingSol(x,y,U,compBox,prmtrs,bt1,bt2,al,c,theta(end),derivative,nh);
prmtrs.h = nh;
prmtrs.tau = tauVector(end);
save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
    '_h0' num2str(nh * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);

PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
PrintResults(solutionNorms,mu);
PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c);
return;
% Enlarge Domain from existing save:
%old domain parameters from the save must be the same as define here:
x_st = -40.0;    y_st = -40.0;
x_end = 40.0;    y_end = 40.0;
%new domain boundaries:
x_st2 = -80.0;   y_st2 = -80.0;
x_end2 = 80.0;   y_end2 = 80.0;

compBox = struct('x_st',{x_st},'x_end',{x_end},'y_st',{y_st},'y_end',...
    {y_end},'x_st2',{x_st2},'x_end2',{x_end2},'y_st2',{y_st2},'y_end2',{y_end2});
tic
if(length(tauVector)<iterMax && UseExtendedDomain == 1 && size(bigUTimeDerivative,1)~=1)
    fprintf('\nLarge Domamin Calculations:\n\n');
    prmtrs.checkBoundary = 0;
    prmtrs.eps = 1.0e-14;
    prmtrs.plotResidual = 0;
    prmtrs.tau = tauVector(end);
    [bigU,bigUTimeDerivative,P,U,newBigIC,solutionNorms,theta,c1,c2,tauVector,angl,sw_div] =...
    PrepareICForEnlargedDomain(bigU,compBox,prmtrs,al,bt1,bt2,c,c1,theta(end),derivative);
    x=x_st2:h:x_end2; y=y_st2:h:y_end2;
end
if(sw_div == 1)
    return;
end
save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
  '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);

fprintf('elapsed time = %d \n', toc);
PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
PrintResults(solutionNorms,mu);
PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c);
return;

% Create two waves from one existing wave:
lastTheta=theta(end); 
last_tau = tauVector(end); 

LeftHalfU = U(:,2:401);
lastU = [fliplr(LeftHalfU) U(:,1:401)] + U;
clear('LeftHalfU');
LeftHalfP = P(:,2:401);
lastP = [fliplr(LeftHalfP) P(:,1:401)] + P;
clear('LeftHalfP');

[zeroX,zeroY]=GetZeroNodes(x,y);
xx=x(zeroX:end); yy=y(zeroY:end); 
    
figure(2);
mesh(xx,yy,( U'))%bbU'
xlabel('x');    ylabel('y');
title('lastU')
    
[bigU,bigUTimeDerivative,P,U,thetaCont,mu,solutionNormsCont,tauVecCont,anglCont,sw_div] =...
       sol_ch_v8(lastU,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastP);
if(sw_div == 1)
    return;
end

% Continue from last iteration:
lastTheta=theta(end); lastU=U; lastP = P;  last_tau = tauVector(end); 
if(prmtrs.useZeroBoundary == 2)
    [bigU,bigUTimeDerivative,P,U,thetaCont,mu,solutionNormsCont,tauVecCont,anglCont,sw_div] =...
        sol_ch_v8ZeroBnd(lastU,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastP);
else
    [bigU,bigUTimeDerivative,P,U,thetaCont,mu,solutionNormsCont,tauVecCont,anglCont,sw_div] =...
        sol_ch_v8(lastU,x,y,prmtrs,bt1,bt2,al,c,lastTheta,derivative,lastP);
end
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
solutionNorms.UL2Norm = solutionNormsCont.UL2Norm;
solutionNorms.PL2Norm = solutionNormsCont.PL2Norm;
solutionNorms.boundaryFunctionUvsUInfNorm = solutionNormsCont.boundaryFunctionUvsUInfNorm;
solutionNorms.BoundaryFunctionPvsPInfNorm = solutionNormsCont.BoundaryFunctionPvsPInfNorm;
solutionNorms.UInfNorm = solutionNormsCont.UInfNorm;
solutionNorms.PInfNorm = solutionNormsCont.PInfNorm;

save (['SavedWorkspaces\' GetICName(ICSwitch) 'IC_2W_' num2str(floor(x_end2)) '_ZB'  num2str(useZeroBoundary) '_bt' num2str(bt) '_c0' num2str(floor(c*100)) ...
    '_h0' num2str(h * 100,'%.02d') '_O(h^' num2str(  size( secondDerivative, 2 ) - 1  ) ')']);

PlotResidualInfNormTauAndUvsUpInfNorm(solutionNorms,tauVector,angl);
PrintResults(solutionNorms,mu);
PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c);
return;

DrawSolution(x,y,h,al,bt,c,theta,bigU,bigUTimeDerivative,bigIC,U,compBox,secondDerivative);

DrawAsymptoticTerms(x,y,h,al,bt,c,theta,bigU,bigUTimeDerivative,bigIC,U,compBox,secondDerivative);

PlotAssymptVsSolu( x, y, h, bigU, mu.muU*theta(end), c/sqrt(bt), 2 );
PlotAssymptotics(x,y,h,zeroX,zeroY,bigU);
DrawDerivativesOfSolution(bigU,compBox,x,y,h,zeroX,zeroY,c,derivative);
return;



    [X,Y]=Domain(x,y);

    c12 = 1-c^2;
    %newBigIC = c1*lastTheta*(c12*X.^2-Y.^2)./(c12*X.^2+Y.^2).^2;
    newBigIC =  X.*(c12*X.^2-3*Y.^2)./(c12*X.^2+Y.^2).^3; 
    %newBigIC =  (c12^2* X.^4 - 6*c12 * X.^2 .* Y.^2 + Y.^4)./(c12*X.^2+Y.^2).^4;
    
    fig_ss13=figure(13);
    set(fig_ss13, 'OuterPosition', [0.0      	30.0        380.0     340.0]);
    %*sqrt(1-c^2)
    mesh(x,y,newBigIC');
    xlabel('x');    ylabel('y');
    title('end solution')
    %axis([x_st2 x_end2 y_st2 y_end2 -0.5 1]);
    axis([x_st x_end y_st y_end -0.01 .01]);
    colorbar;
    caxis([-0.001 .001]);
    view(0,90);
    


        hx = h;
        hy = h;
        result = hx*hy*sum( sum( bigU(2:end-1, 2:end-1 ) ) ) +...
            hx*hy/2 * ( sum( bigU(1, 2:end-1 ) ) + sum( bigU(end, 2:end-1 ) ) )+...
            hx*hy/2 * ( sum( bigU(2:end-1, 1 ) ) + sum( bigU(2:end-1, end ) ) )+...
            hx*hy/4 * ( bigU(1, 1) + bigU(1, end) + bigU(end, 1) + bigU(end, end) )

figure(1)
mesh(x(zeroX:end),y(zeroY:end),UTilda');
figure(2)
mesh(x(zeroX:end),y(zeroY:end),Uup');
figure(3)
mesh(x(zeroX:end),y(zeroY:end),Uup1');        
