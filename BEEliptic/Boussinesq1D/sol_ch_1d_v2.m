function [bigU,bigUTimeDerivative,PsiUp,PhiUp,thetaVector,mu,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_1d_v2(Phi,x,y,prmtrs,bt1,bt2,al,c,theta,derivative,Psi)

    sw = 0;  
    if (nargin == 11) 
        sw=1; 
    end

    % constants
    tau = prmtrs.tau;
    h=prmtrs.h;
    eps = prmtrs.eps;
    tauMax = 10*tau;
    tauDecreasedIteration = 1;
    iterMax = prmtrs.iterMax;
    bt = bt1/bt2;
    autoStop = 0;
    afterCounter = 5000;

    manualStop=0;
    boundaryHit = 0;
    sw_div = 0;

    UvsUupInfNorm = ones(1,iterMax);
    thetaVector = zeros(1,iterMax);
    tauVector = zeros(1,iterMax); 
    tauVector(1) = tau;
    residualInfNorm = ones(1,iterMax/10);
    angl = zeros(1,iterMax/10);
    step = Step(h);
    zeroMatrix = zeros(size(Phi));

    Phix = XDerivativeEvenFunctions2D(Phi, zeroMatrix, derivative.first);
    Phixy = YDerivativeEvenFunctions2D(Phix, zeroMatrix, derivative.first);
    deltaPhi = DeltaEvenFunctions2D(Phi, zeroMatrix,derivative.second);
    if(sw == 0)
        Psi = (bt*c^2-1)*deltaPhi - 2*(bt*c^2+1)*Phixy + theta*al*bt*Phi.^2;        
    end
     
	iterCounter=1;
	subCounter=1;
	PhiSquare = Phi .^2; 
    % check that 2*(bt*c^2+1)*Phixy(1,1) == 0
	thetaVector(iterCounter) = ( (bt*c^2-1)*deltaPhi(1,1) - 2*(bt*c^2+1)*Phixy(1,1) - Psi(1,1) )/(al*bt*PhiSquare(1,1));
	residualInfNorm(1) =  GetResidualInfNorm(...
         al,...
         bt,...
         c,...
         thetaVector,...
         iterCounter,...
         Phi,...
         PhiSquare,...
         Phixx+Phiyy,...
         zeroMatrix,...
         muU*outerRigthBoundaryF,...
         muU*outerTopBoundaryF,...
         derivative.second);
     
	PsiUp = Psi;
	PhiUp = Phi;
    minResidual = min(1000,residualInfNorm(subCounter));
    crrntResidual = 0;
    fig9=figure(9);
    set(gcf,'currentcharacter','T');
    flag = 0; ax = 1; ay = 1; 
    while( 1 > 0 ) %~( residualInfNorm(subCounter) < eps )
        Psi = PsiUp;
        Phi = PhiUp;
        iterCounter=iterCounter+1;

        %deltaU = DeltaEvenFunctions(U, zeroMatrix, muU*outerRigthBoundaryF, muU*outerTopBoundaryF,derivative.second); %<--- TIME CONSummable
        PhiSquare = Phi .^2;
        Phixx = XDerivativeEvenFunctions2D(Phi, zeroMatrix,muU*outerRigthBoundaryF,derivative.second);
        Phiyy = YDerivativeEvenFunctions2D(Phi, zeroMatrix,muU*outerTopBoundaryF,derivative.second);
         % TIME CONSumable       |
         %                      \|/
         %                       '
         
        thetaVector(iterCounter) = ( (bt*c^2-1)*( Phixx(1,1)+Phiyy(1,1) ) - Psi(1,1) )/(al*bt*PhiSquare(1,1));
        PsiUp = Psi + (tau)*( DeltaEvenFunctions2D(Psi - bt*(c^2-1)*Phi, zeroMatrix,derivative.second) +...
            XDerivativeEvenFunctions2D( YDerivativeEvenFunctions2D( 2*Psi + 2*bt*(c^2+1)*Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first) );  

        PhiUp = Phi + tau*( (bt*c^2 - 1) * (Phixx+Phiyy) - PsiUp + al*bt*thetaVector(iterCounter) * PhiSquare +...
             - 2*(bt*c^2 + 1) * XDerivativeEvenFunctions2D( YDerivativeEvenFunctions2D( Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first) );
        
        UvsUupInfNorm(iterCounter) = max(max(abs(Phi-PhiUp)));
        
        if(mod(iterCounter,10) ==0)        
           
           crrntResidual = Get2DResidual(al,bt,c,thetaVector,iterCounter,Phi,PhiSquare,Phixx+Phiyy,...
               zeroMatrix,muU*outerRigthBoundaryF,muU*outerTopBoundaryF,derivative.second);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=abs(thetaVector(iterCounter))*max(max(abs(crrntResidual(1:end-8,1:end-8))));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           [flag, Px, Py] = StopCriteria(x, y, zeroX, zeroY, Phi, ax, ay, minResidual, eps);
           
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               %fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n eps = %.15e\n', abs( ax-Px(1) ), abs( ay-Py(1) ), eps );
                
               if(prmtrs.plotResidual)
                   PlotResidual(x(zeroX:end),y(zeroY:end),crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x,y,zeroX, Phi, muU*outerTopBoundaryF);
               end
               if(prmtrs.plotAssympt)
                   PlotAssymptVsSolu( x(zeroX:end), y(zeroY:end), h, 1, 1, Phi, muU*thetaVector(iterCounter), c);
               end
               
               manualStop = Stop(fig9,manualStop);
           end
           ax = Px(1);
           ay = Py(1);
           boundaryHit = IsBoundaryHit(checkBnd,crrntResidual*thetaVector(iterCounter),...
               residualInfNorm,subCounter,derivative.second);
           sw_div = IsAlgoDiverging(subCounter,residualInfNorm); 
           
           angl(subCounter) =  Deviation(residualInfNorm,subCounter);
        end
        
        tauVector(iterCounter) = tau;        
        [tau, tauMax, tauIncreasedIteration, tauDecreasedIteration] =...
        DefineCurrentTau(prmtrs.tau, subCounter, iterCounter, iterMax, tau,  tauMax, tauIncreasedIteration,...
           tauDecreasedIteration, residualInfNorm,UvsUupInfNorm, minResidual, angl , manualStop, flag, crrntResidual);
        
        if( flag == 1 )
           afterCounter = afterCounter - 1; 
        end
        if(afterCounter <= 0)
            autoStop = 1;
        end
        if(sw_div ==1 || manualStop ==2 || boundaryHit ==1 || iterCounter ==iterMax || autoStop ==1)
            break;
        end
    end
    UvsUupInfNorm = thetaVector(iterCounter)*UvsUupInfNorm(1:iterCounter);
    tauVector = tauVector(1:iterCounter);
    thetaVector = thetaVector(1:iterCounter);
    angl = angl(1:subCounter);
    [muU, muP]=FindBoundaryConstants(PhiUp,PsiUp,innerBoundaryUF,innerBoundaryPF,step);
    mu = struct('muU',{muU},'muP',{muP});

    [bigU,bigUTimeDerivative] = GetBigSolution(x,y,c,bt,thetaVector(iterCounter),...
        PhiUp,bigZeroMatrix,derivative.first,outerTopBoundaryF*muU);
    
    [solutionNorms] = CalculateSolutionNorms(Phi,PhiUp,Psi,PsiUp,UvsUupInfNorm,crrntResidual,...
        residualInfNorm,muU*innerBoundaryUF,muP*innerBoundaryPF,subCounter,thetaVector(iterCounter),step,h);
end

%z1 = YDerivativeEvenFunctions2(P,zeroMatrix,yDerBnd,derivative.second);
%z2 = YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second);
%mesh(x(zeroX:end),y(zeroY:end),(z2-z1)');
%DrawResidual(x(zeroX:end),y(zeroY:end),al,bt,c,thetaVector(iterCounter),c1,U,derivative.second);
