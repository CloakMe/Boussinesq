function [PhiUp,PsiUp,thetaVector,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_1d_v3(Phi,x,y,prmtrs,bt1,bt2,al,c,theta,derivative,Psi)

    %equation system after VC s = x-ct, r = x+ct and removing u_xxtt (btExt = 0) transforms to
    % pSi = (btExt*bt*c^2-1) * Delta(phi) - 2*(btExt*bt*c^2+1) * phi_sr + al*bt*phi^2
    % bt(c^2-1) * Delta(phi) - 2*bt*(c^2+1) * phi_sr = Delta(pSi) + 2 * pSi_sr
    % which transforms to (btExt = 0):
    % pSi = (-1) * Delta(phi) - 2*(+1) * phi_sr + al*bt*phi^2
    % bt(c^2-1) * Delta(phi) - 2*bt*(c^2+1) * phi_sr = Delta(pSi) + 2 * pSi_sr    
    
    sw = 0;  
    if (nargin == 11) 
        sw=1; 
    end

    % constants
    btExt = 0;
    tau = prmtrs.tau;
    h=prmtrs.h;
    eps = prmtrs.eps;
    tauMax = 10*tau;
    tauIncreasedIteration = 60;
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

    Phix = XDerivative(Phi, zeroMatrix, derivative.first);
    Phixy = YDerivative(Phix, zeroMatrix, derivative.first);
    deltaPhi = XDerivative(Phi, zeroMatrix, derivative.second) + YDerivative(Phi, zeroMatrix, derivative.second);
    if(sw == 0)
        Psi = - deltaPhi - 2*Phixy + theta*al*bt*Phi.^2;        
    end
     
	iterCounter=1;
	subCounter=1;
	PhiSquare = Phi .^2;
    % check that 2*(bt*c^2+1)*Phixy(1,1) == 0
	thetaVector(iterCounter) = ( Psi(zeroX,zeroY) - (-1)*deltaPhi(zeroX,zeroY) + 2*(+1)*Phixy(zeroX,zeroY) )/(al*bt*PhiSquare(zeroX,zeroY));
	residual = GetResidual(prmtrs.type, bt, btExt, c, Phi, al*bt*thetaVector(iterCounter)*PhiSquare, deltaPhi, zeroMatrix, derivative);
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
    %crossSect = max(1,ceil(length(x)/8));
    while( 1 > 0 ) %~( residualInfNorm(subCounter) < eps )
        iterCounter=iterCounter+1;

        PhiSquare = Phi .^2;
        deltaPhi = XDerivative(Phi, zeroMatrix, derivative.second) + YDerivative(Phi, zeroMatrix, derivative.second);
         % TIME CONSumable       |
         %                      \|/
         %                       '         
        thetaVector(iterCounter) = ( Psi(zeroX,zeroY) - (-1)*deltaPhi(zeroX,zeroY) + 2*(+1)*Phixy(zeroX,zeroY) )/(al*bt*PhiSquare(zeroX,zeroY));
        PsiUp = Psi + tau * ( XDerivative(Psi - bt*(c^2-1)*Phi, zeroMatrix,derivative.second) + YDerivative(Psi - bt*(c^2-1)*Phi, zeroMatrix,derivative.second) +...
            XDerivative( YDerivative( 2*Psi + 2*bt*(c^2+1)*Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first) );  

        PhiUp = Phi + tau * ( (- 1) * deltaPhi - PsiUp + al*bt*thetaVector(iterCounter) * PhiSquare -...
            2*(+ 1) * XDerivative( YDerivative( Phi, zeroMatrix, derivative.first), zeroMatrix, derivative.first) );
                
        UvsUupInfNorm(iterCounter) = max(max(abs(Phi-PhiUp)));
        
        if(mod(iterCounter,10) ==0)        
           
           crrntResidual = GetResidual(prtmrs.type, bt, btExt, c, Phi, al*bt*thetaVector(iterCounter)*PhiSquare, deltaPhi, zeroMatrix, derivative);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=abs(thetaVector(iterCounter))*max(max(abs(crrntResidual)));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           %[flag, Px, Py] = StopCriteria(x, y, zeroX, zeroY, Phi, ax, ay, minResidual, eps);
           
%            ratioPhi = abs(Phi(crossSect,zeroY))/abs(Phi(crossSect*2,zeroY));
%            if(ratioPhi > 1)
%                ratioPhi = 1/ratioPhi;
%                warning('ratioPhi > 1 ');
%            end
%            Phi(1,:) = ratioPhi*Phi(crossSect,:);
%            Phi(end,:) = ratioPhi*Phi(end-crossSect,:);
%            Phi(:,1) = ratioPhi*Phi(:,crossSect);
%            Phi(:,end) = ratioPhi*Phi(:,end-crossSect);
%            
%            ratioPsi = abs(Psi(crossSect*2,zeroY))/abs(Psi(crossSect*4,zeroY));
%            if(ratioPsi > 1)
%                ratioPsi = 1;
%                warning('ratioPsi > 1');
%            end
%            Psi(1,:) = ratioPsi*Psi(crossSect*2,:);
%            Psi(end,:) = ratioPsi*Psi(end-crossSect*2,:);
%            Psi(:,1) = ratioPsi*Psi(:,crossSect*2);
%            Psi(:,end) = ratioPsi*Psi(:,end-crossSect*2);
           
%            figure(20)
%            plot(x,Psi(1,:),'b', x, Psi(crossSect,:), 'r', x, Psi(crossSect*2,:), 'k');
%            xlabel('x');    ylabel('psi');
%            title('Boundary (Psi)');
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               %fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n eps = %.15e\n', abs( ax-Px(1) ), abs( ay-Py(1) ), eps );
                
               if(prmtrs.plotResidual)
                   PlotResidual(x, y, crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x, y, 1, Phi);
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
        [tau, tauMax, tauIncreasedIteration, tauDecreasedIteration] =...
            DefineCurrentTau(prmtrs.tau, subCounter, iterCounter, iterMax, tau,  tauMax, tauIncreasedIteration,...
                tauDecreasedIteration, residualInfNorm,UvsUupInfNorm, minResidual, angl , manualStop, flag, crrntResidual);
        
        if( flag == 1 )
           afterCounter = afterCounter - 1; 
        end
        if(afterCounter <= 0)
            autoStop = 1;
        end
        if(sw_div == 1 || manualStop == 2 || iterCounter == iterMax || autoStop == 1)
            break;
        end
        
        Psi = PsiUp;
        Phi = PhiUp;        
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
