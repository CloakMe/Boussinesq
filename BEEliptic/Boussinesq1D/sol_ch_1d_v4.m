function [UUp,PUp,thetaVector,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_1d_v4(U,x,t,prmtrs,bt1,bt2,al,c,theta,derivative,P)

    sw = 0;  
    if (nargin == 11) 
        sw=1; 
    end

    % constants
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
    [zeroX,zeroT]=GetZeroNodes(x,t);
    
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
    zeroMatrix = zeros(size(U));
    
    Uxx = XDerivative(U, zeroMatrix, derivative.second);
    hypU = YDerivative(U, zeroMatrix, derivative.second) - Uxx;
    if(sw == 0)
        P = theta*al*bt*U.^2 - Uxx;        
    end
     
	iterCounter=1;
	subCounter=1;
	USquare = U .^2;
    % check that 2*(bt*c^2+1)*Phixy(1,1) == 0
	thetaVector(iterCounter) = ( P(zeroX,zeroT) + Uxx(zeroX,zeroT))/(al*bt*USquare(zeroX,zeroT));
	residual = GetResidual(prmtrs.type, bt, c, U, al*thetaVector(iterCounter), hypU, zeroMatrix, derivative);
%     figure(11)
%     mesh(x,t,residual*thetaVector(iterCounter));
%     title('residual');
    residualInfNorm(1) = abs( thetaVector(iterCounter) ) * max(max(abs(residual)));
	PUp = P;
	UUp = U;
    minResidual = min(1000,residualInfNorm(subCounter));
    crrntResidual = 0;
    fig9=figure(9);
    set(gcf,'currentcharacter','T');
    flag = 0; ax = 1; ay = 1; 
    %crossSect = max(1,ceil(length(x)/8));
    while( 1 > 0 ) %~( residualInfNorm(subCounter) < eps )
        iterCounter=iterCounter+1;

        USquare = U .^2;        
        Uxx = XDerivative(U, zeroMatrix, derivative.second);
         % TIME CONSumable       |
         %                      \|/
         %                       '         
        thetaVector(iterCounter) = ( P(zeroX,zeroT) + Uxx(zeroX,zeroT))/(al*bt*USquare(zeroX,zeroT));
        PUp = P + tau * (bt * YDerivative(U, zeroMatrix, derivative.second) - bt * XDerivative(U, zeroMatrix,derivative.second) -...
            XDerivative( P, zeroMatrix, derivative.first) );  

        UUp = U + tau * (PUp + XDerivative(U, zeroMatrix, derivative.second) - al*bt*thetaVector(iterCounter) * U .^2);
                
        UvsUupInfNorm(iterCounter) = max(max(abs(U-UUp)));
        
        if(mod(iterCounter,10) ==0)        
           
           crrntResidual = GetResidual(prmtrs.type, bt, c, U, al*thetaVector(iterCounter), 0, zeroMatrix, derivative);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=abs(thetaVector(iterCounter))*max(max(abs(crrntResidual)));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           %[flag, Px, Py] = StopCriteria(x, y, zeroX, zeroY, U, ax, ay, minResidual, eps);
           
%          figure(20)
%          plot(x,Psi(1,:),'b', x, Psi(crossSect,:), 'r', x, U(crossSect*2,:), 'k');
%          xlabel('x');    ylabel('psi');
%          title('Boundary (Psi)');
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               %fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n eps = %.15e\n', abs( ax-Px(1) ), abs( ay-Py(1) ), eps );
                
               if(prmtrs.plotResidual)
                   PlotResidual(x, t, crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x, t, 1, U);
               end
               if(prmtrs.plotAssympt)
                   PlotAssymptVsSolu( x(zeroX:end), t(zeroT:end), h, 1, 1, U, thetaVector(iterCounter), c);
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
        
        P = PUp;
        U = UUp;        
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
