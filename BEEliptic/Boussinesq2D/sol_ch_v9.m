function [bigU,bigUTimeDerivative,Pup,Uup,thetaVector,c1,c2,solutionNorms,tauVector, angl]=...
    sol_ch_v9(U,x,y,prmtrs,bt1,bt2,al,c,theta,zeroX,zeroY,derivative,P)

sw = 0;  
if (nargin == 13) 
    sw=1; 
end

% constants
   tau = prmtrs.tau;
   h=prmtrs.h;
   eps = prmtrs.eps;
   checkBnd = prmtrs.checkBoundary;
   tauMax = 10*tau;
   tauIncreasedIteration = 60;
   tauDecreasedIteration = 1;
   iterMax = prmtrs.iterMax;
   bt = bt1/bt2;
   
   sx = length(x);
   sy = length(y);
   stopSwitch=0;
   boundaryHit = 0;
   sw_div = 0;

   UvsUupInfNorm = ones(1,iterMax);
   thetaVector = zeros(1,iterMax);
   tauVector = zeros(1,iterMax); 
   tauVector(1) = tau;
   residualInfNorm = ones(1,iterMax/10);
   angl = zeros(1,iterMax/10);
   step = Step(h);
   
   % approxBoundaryF is over the augmented domain 
   % where four points were added top and bottom
   approxBoundaryF = GetApproximationForBoundary(x(zeroX:end),y(zeroY:end),h,c/sqrt(bt));
   
   outerTopBoundaryF=approxBoundaryF(1:end-4,end-3:end); 
   outerRigthBoundaryF=approxBoundaryF(end-3:end,1:end-4); 
   
   xLength=length(x(zeroX:end));  yLength=length(y(zeroY:end));
 
   [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

   innerBoundaryUF=[approxBoundaryF(xHeiA,yLen)  approxBoundaryF(xLen,yHeiA)'];
   innerBoundaryPF = (bt-c^2)*innerBoundaryUF;
   
   bigZeroMatrix = zeros(sx,sy);
   zeroMatrix = zeros(size(U));
    c1 = FindBoundaryConstants(U,0*U,innerBoundaryUF,innerBoundaryPF,step);
    deltaU = DeltaEvenFunctions(U,zeroMatrix,c1*outerRigthBoundaryF,c1*outerTopBoundaryF,derivative.second);
    if(sw == 0)
        P = (bt-c^2)* U - (1 - c^2)*deltaU + theta*al*bt*U.^2;
    end
     
     iterCounter=1;
     subCounter=1;
     Usquare = U .^2;
     thetaVector(iterCounter) = (P(1,1) -  (bt-c^2)*U(1,1) + (1 - c^2) * deltaU(1,1) )/(al*bt*Usquare(1,1));
     residualInfNorm(1) =  GetResidualInfNorm(al,bt,c/sqrt(bt),thetaVector,0,iterCounter,U,Usquare,...
         deltaU,zeroMatrix,c1*outerRigthBoundaryF,c1*outerTopBoundaryF,derivative.second);
     
     Pup = P;
     Uup = U;
     minResidual = min(1000,residualInfNorm(subCounter));
     crrntResidual = 0;
     fig9=figure(9);
     set(gcf,'currentcharacter','T');
     flag = 0; ax = 1; ay = 1; 
    while( 1 > 0 ) %~( residualInfNorm(subCounter) < eps )
        P = Pup;
        U = Uup;
        iterCounter=iterCounter+1;
        [c1,c2] = FindBoundaryConstants(U,P,innerBoundaryUF,innerBoundaryPF,step);
        deltaU = DeltaEvenFunctions(U,zeroMatrix,c1*outerRigthBoundaryF,c1*outerTopBoundaryF,derivative.second); %<--- TIME CONSummable
        Usquare = U .^2;
        thetaVector(iterCounter) = (P(1,1) -  (bt-c^2)*U(1,1) + (1 - c^2) * deltaU(1,1) )/(al*bt*Usquare(1,1));
         % TIME CONSumable       |
         %                      \|/
         %                       '
        yDerBnd = c2*(bt-c^2)*outerTopBoundaryF;
        xDerBnd = (c^2*c1 + c2*(bt-c^2))*outerRigthBoundaryF;
        
        Pup = P + (tau)*(YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second) +...
            XDerivativeEvenFunctions(P-c^2*(deltaU - U),zeroMatrix,xDerBnd,derivative.second));  
        
        Uup = U + tau*( Pup - (al*bt*thetaVector(iterCounter) )*Usquare +...
            ((1 - c^2))* deltaU - (bt - c^2)*U );
        
        if(mod(iterCounter,6000) ==0)

        end
        UvsUupInfNorm(iterCounter) = max(max(abs(U-Uup)));
        
        if(mod(iterCounter,10) ==0)
           
           [flag, axNew, ayNew] = StopCriteria(x, y, h, zeroX, zeroY, U, ax, ay, eps);
           crrntResidual = Get2DResidual(al,bt,c/sqrt(bt),thetaVector,c1,iterCounter,U,Usquare,deltaU,...
               zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,derivative.second);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=thetaVector(iterCounter)*max(max(abs(crrntResidual(1:end-8,1:end-8))));
           angl(subCounter) =  Deviation(residualInfNorm,subCounter);
           minResidual = min(minResidual,residualInfNorm(subCounter));
           if(minResidual > .001)
               flag = 0;
           end
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n', abs( ax-axNew ), abs( ay-ayNew ) );
               
               if(prmtrs.plotResidual)
                   PlotResidual(x(zeroX:end),y(zeroY:end),crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x,y,zeroX, U, outerTopBoundaryF,c1);
               end
               if(prmtrs.plotAssympt)
                   PlotAssymptotics(x,y,h,zeroX,zeroY,U);
               end
               
               stopSwitch = Stop(fig9,stopSwitch);
           end
           ax = axNew;
           ay = ayNew;
           boundaryHit = IsBoundaryHit(checkBnd,crrntResidual*thetaVector(iterCounter),...
               residualInfNorm,subCounter,derivative.second);
           sw_div = IsAlgoDiverging(subCounter,residualInfNorm); 
        end
           
        tauVector(iterCounter) = tau;
        
        [tau, tauMax, tauIncreasedIteration, tauDecreasedIteration] =...
        DefineCurrentTau(subCounter, iterCounter, iterMax, tau,  tauMax, tauIncreasedIteration,...
            tauDecreasedIteration, residualInfNorm,UvsUupInfNorm, minResidual, angl , stopSwitch,crrntResidual);
    
        if(sw_div ==1 || stopSwitch ==2 || boundaryHit ==1 || iterCounter ==iterMax || flag ==1)
            break;
        end
    end
    UvsUupInfNorm = thetaVector(iterCounter)*UvsUupInfNorm(1:iterCounter);
    tauVector = tauVector(1:iterCounter);
    thetaVector = thetaVector(1:iterCounter);
    angl = angl(1:subCounter);
    [c1,c2]=FindBoundaryConstants(Uup,Pup,innerBoundaryUF,innerBoundaryPF,step);
    
    [bigU,bigUTimeDerivative] = GetBigSolution(x,y,c,bt,c1,thetaVector(iterCounter),...
        Uup,bigZeroMatrix,derivative.first,outerTopBoundaryF);
   
    [solutionNorms] = CalculateSolutionNorms(U,Uup,P,Pup,UvsUupInfNorm,crrntResidual,...
        residualInfNorm,innerBoundaryUF,innerBoundaryPF,subCounter,thetaVector(iterCounter),step,h,c1,c2);
end

%z1 = YDerivativeEvenFunctions2(P,zeroMatrix,yDerBnd,derivative.second);
%z2 = YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second);
%mesh(x(zeroX:end),y(zeroY:end),(z2-z1)');
%DrawResidual(x(zeroX:end),y(zeroY:end),al,bt,c,thetaVector(iterCounter),c1,U,derivative.second);
