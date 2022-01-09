function [PhiUp,PsiUp,thetaVector,solutionNorms,tauVector,angl,sw_div]=...
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
    [zeroX,zeroY]=GetZeroNodes(x,y);
    
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
	thetaVector(iterCounter) = ( Psi(zeroX,zeroY) - (bt*c^2-1)*deltaPhi(zeroX,zeroY) + 2*(bt*c^2+1)*Phixy(zeroX,zeroY) )/(al*bt*PhiSquare(zeroX,zeroY));
	residual = GetResidual(bt, c, Phi, al*bt*thetaVector(iterCounter)*PhiSquare, deltaPhi, zeroMatrix, derivative);
    figure(11)
    mesh(x,y,residual*thetaVector(iterCounter));
    title('residual');
    residualInfNorm(1) = abs( thetaVector(iterCounter) ) * max(max(abs(residual)));
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

        PhiSquare = Phi .^2;
        deltaPhi = DeltaEvenFunctions2D(Phi, zeroMatrix,derivative.second);
         % TIME CONSumable       |
         %                      \|/
         %                       '
         
        thetaVector(iterCounter) = ( Psi(zeroX,zeroY) - (bt*c^2-1)*deltaPhi(zeroX,zeroY) + 2*(bt*c^2+1)*Phixy(zeroX,zeroY) )/(al*bt*PhiSquare(zeroX,zeroY));
        PsiUp = Psi + tau * ( DeltaEvenFunctions2D(Psi - bt*(c^2-1)*Phi, zeroMatrix,derivative.second) +...
            XDerivativeEvenFunctions2D( YDerivativeEvenFunctions2D( 2*Psi + 2*bt*(c^2+1)*Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first) );  

        PhiUp = Phi + tau * ( (bt*c^2 - 1) * deltaPhi - PsiUp + al*bt*thetaVector(iterCounter) * PhiSquare +...
             - 2*(bt*c^2 + 1) * XDerivativeEvenFunctions2D( YDerivativeEvenFunctions2D( Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first) );
         
        PsiUp(1,:) = Psi(1,:);
        PsiUp(end,:) = Psi(end,:);
        PsiUp(:,1) = Psi(:,1);
        PsiUp(:,end) = Psi(:,end);
        PhiUp(1,:) = Phi(1,:);
        PhiUp(end,:) = Phi(end,:);
        PhiUp(:,1) = Phi(:,1);
        PhiUp(:,end) = Phi(:,end);
        
        UvsUupInfNorm(iterCounter) = max(max(abs(Phi-PhiUp)));
        
        if(mod(iterCounter,10) ==0)        
           
           crrntResidual = GetResidual(bt, c, Phi, al*bt*thetaVector(iterCounter)*PhiSquare, deltaPhi, zeroMatrix, derivative);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=abs(thetaVector(iterCounter))*max(max(abs(crrntResidual(1:end-8,1:end-8))));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           %[flag, Px, Py] = StopCriteria(x, y, zeroX, zeroY, Phi, ax, ay, minResidual, eps);
           
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               %fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n eps = %.15e\n', abs( ax-Px(1) ), abs( ay-Py(1) ), eps );
                
               if(prmtrs.plotResidual)
                   PlotResidual(x, y, crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x,y,1, Phi, 0);
               end
               if(prmtrs.plotAssympt)
                   PlotAssymptVsSolu( x(zeroX:end), y(zeroY:end), h, 1, 1, Phi, thetaVector(iterCounter), c);
               end
               
               manualStop = Stop(fig9,manualStop);
           end
           sw_div = IsAlgoDiverging(subCounter,residualInfNorm); 
           
           angl(subCounter) =  Deviation(residualInfNorm,subCounter);
        end
        
        tauVector(iterCounter) = tau;        
        %[tau, tauMax, tauIncreasedIteration, tauDecreasedIteration] =...
        %    DefineCurrentTau(prmtrs.tau, subCounter, iterCounter, iterMax, tau,  tauMax, tauIncreasedIteration,...
        %        tauDecreasedIteration, residualInfNorm,UvsUupInfNorm, minResidual, angl , manualStop, flag, crrntResidual);
        
        if( flag == 1 )
           afterCounter = afterCounter - 1; 
        end
        if(afterCounter <= 0)
            autoStop = 1;
        end
        if(sw_div == 1 || manualStop == 2 || iterCounter == iterMax || autoStop == 1)
            break;
        end
    end
    UvsUupInfNorm = thetaVector(iterCounter)*UvsUupInfNorm(1:iterCounter);
    tauVector = tauVector(1:iterCounter);
    thetaVector = thetaVector(1:iterCounter);
    angl = angl(1:subCounter);

    solutionNorms = struct('residualInfNorm',{residualInfNorm}, 'UvsUpInfNorm',{UvsUupInfNorm});
end

%z1 = YDerivativeEvenFunctions2(P,zeroMatrix,yDerBnd,derivative.second);
%z2 = YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second);
%mesh(x(zeroX:end),y(zeroY:end),(z2-z1)');
%DrawResidual(x(zeroX:end),y(zeroY:end),al,bt,c,thetaVector(iterCounter),c1,U,derivative.second);
