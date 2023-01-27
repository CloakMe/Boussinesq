function [bigU,bigUTimeDerivative,Pup,Uup,thetaVector,mu,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_v8ZeroBnd(U,x,y,prmtrs,bt1,bt2,al,c,theta,derivative,P)

    sw = 0;  
    if (nargin == 13) 
        sw=1; 
    end
    
    domainUtils = BEDomainUtils( x, y, length(derivative.first)-1 );
    
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
    autoStop = 0;
    afterCounter = 5000;

    sx = length(x);
    sy = length(y);
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

    % approxBoundaryF is over the augmented domain 
    % where four points were added top and right
    [zeroX,zeroY]=GetZeroNodes(x,y);
    approxBoundaryF = GetApproximationForBoundary(x(zeroX:end),y(zeroY:end),h,c);

    outerTopBoundaryF=approxBoundaryF(1:end-4,end-3:end); 
    outerRigthBoundaryF=approxBoundaryF(end-3:end,1:end-4); 

    xLength=length(x(zeroX:end));  yLength=length(y(zeroY:end));

    [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

    innerBoundaryUF=[approxBoundaryF(xHeiA,yLen)  approxBoundaryF(xLen,yHeiA)'];
    innerBoundaryPF = bt*(1-c^2)*innerBoundaryUF;

    %bigZeroMatrix = zeros(sx,sy);
    zeroMatrix = zeros(size(U));
    if(prmtrs.useZeroBoundary > 0)
            muU = 0; muP = 0;
        else
            [muU,muP] = FindBoundaryConstants(U,P,innerBoundaryUF,innerBoundaryPF,step);
        end
    %deltaU = DeltaEvenFunctions(U, zeroMatrix, muU*outerRigthBoundaryF, muU*outerTopBoundaryF, derivative.second);
    Uxx = domainUtils.XDerivativeEvenFunZeroBnd(U)/h^2;
    Uyy = domainUtils.YDerivativeEvenFunZeroBnd(U)/h^2;
    if(sw == 0)
        %P = bt*(1-c^2)*U - (1-bt*c^2)*(Uyy + Uxx) + theta*al*bt*U.^2;
        P = bt*(1-c^2)*U - (1-bt*c^2)*Uyy - Uxx + theta*al*bt*U.^2;        
    end
     
     iterCounter=1;
     subCounter=1;
     Usquare = U .^2; %
     %thetaVector(iterCounter) = (P(1,1) - bt*(1-c^2)*U(1,1) + (1-bt*c^2)*( Uyy(1,1) + Uxx(1,1)) )/(al*bt*Usquare(1,1));
     thetaVector(iterCounter) = (P(zeroX,zeroY) - bt*(1-c^2)*U(zeroX,zeroY) + (1-bt*c^2)*Uyy(zeroX,zeroY) + Uxx(zeroX,zeroY) )/(al*bt*Usquare(zeroX,zeroY));
     
     mid = floor(length(derivative.second)/2);
     crrntResidual = Get2DResidual(al,bt,c,theta,iterCounter,U,Usquare,...
         Uxx + Uyy,zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,derivative.second, prmtrs, domainUtils );

     residualInfNorm(1)=max(max(abs(crrntResidual(1:end-mid ,1:end-mid ))));
     
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
        if(prmtrs.useZeroBoundary > 0)
            muU = 0; muP = 0;
        else
            [muU,muP] = FindBoundaryConstants(U,P,innerBoundaryUF,innerBoundaryPF,step);
        end
        %deltaU = DeltaEvenFunctions(U, zeroMatrix, muU*outerRigthBoundaryF, muU*outerTopBoundaryF,derivative.second); %<--- TIME CONSummable
        Usquare = U .^2;
        Uxx = domainUtils.XDerivativeEvenFunZeroBnd(U)/h^2;
        Uyy = domainUtils.YDerivativeEvenFunZeroBnd(U)/h^2;
         % TIME CONSumable       |
         %                      \|/
         %                       '
        yDerBnd = muP*bt*(1-c^2)*outerTopBoundaryF;
        xDerBnd = (muU*bt*c^2 + muP*bt*(1-c^2))*outerRigthBoundaryF;
        if( false )
            thetaVector(iterCounter) = (P(1,1) - bt*(1-c^2)*U(1,1) + (1-bt*c^2)*( Uxx(1,1)+Uyy(1,1) ) )/(al*bt*Usquare(1,1));
            Pup = P + (tau)*(domainUtils.YDerivativeEvenFunZeroBnd(P)/h^2 +...
                domainUtils.XDerivativeEvenFunZeroBnd(P-bt*c^2*(Uxx + Uyy - U))/h^2);  
        
            Uup = U + tau*( Pup - (al*bt*thetaVector(iterCounter) )*Usquare +...
                ((1 - bt*c^2))* (Uxx+Uyy) - (bt*(1 - c^2))*U );
        else 
            thetaVector(iterCounter) = (P(zeroX,zeroY) - bt*(1-c^2)*U(zeroX,zeroY) + (1-bt*c^2)*Uyy(zeroX,zeroY) + Uxx(zeroX,zeroY) )/(al*bt*Usquare(zeroX,zeroY));
            Pup = P + (tau)*(domainUtils.YDerivativeEvenFunZeroBnd(P)/h^2 +...
                domainUtils.XDerivativeEvenFunZeroBnd(P)/h^2 + bt*c^2*Uxx );  

            Uup = U + tau*( Pup - (al*bt*thetaVector(iterCounter) )*Usquare +...
                (1 - bt*c^2)*Uyy + Uxx - (bt*(1 - c^2))*U  );        
        end
        
        UvsUupInfNorm(iterCounter) = max(max(abs(U-Uup)));
        figure(2); mesh(x(zeroX:end),y,U');
        if(mod(iterCounter,10) ==0)        
           
           crrntResidual = Get2DResidual(al,bt,c,thetaVector,iterCounter,U,Usquare,Uxx+Uyy,...
               zeroMatrix,muU*outerRigthBoundaryF,muU*outerTopBoundaryF,derivative.second, prmtrs, domainUtils);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=abs(thetaVector(iterCounter))*max(max(abs(crrntResidual(1:end-8,1:end-8))));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           [flag, Px, Py] = StopCriteria(x, y, zeroX, zeroY, U, ax, ay, minResidual, eps);
           
           if(mod(iterCounter,50) ==0)
               figure(2); mesh(x(zeroX:end),y,U');
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               %fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n eps = %.15e\n', abs( ax-Px(1) ), abs( ay-Py(1) ), eps );
               %if(sum(tauVector(1:iterCounter)) > 10455)
               %    autoStop = 1;
               %    fprintf('stop at sum(tauVector) = %.4f \n', sum(tauVector(1:iterCounter)));
               %end
               if(prmtrs.plotResidual)
                   PlotResidual(x(zeroX:end),y(zeroY:end),crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x,y,zeroX, U, muU*outerTopBoundaryF);
               end
               if(prmtrs.plotAssympt)
                   PlotAssymptVsSolu( x(zeroX:end), y(zeroY:end), h, 1, 1, U, muU*thetaVector(iterCounter), c);
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
            figure(2); mesh(x(zeroX:end),y,U');
            break;
        end
    end
    UvsUupInfNorm = thetaVector(iterCounter)*UvsUupInfNorm(1:iterCounter);
    tauVector = tauVector(1:iterCounter);
    thetaVector = thetaVector(1:iterCounter);
    angl = angl(1:subCounter);
    
    
    
	domainUtils1Der = BEDomainUtils( x, y, length(derivative.first)-1, 1 ); 
	if( prmtrs.useZeroBoundary == 2 )
        [muU, muP]=FindBoundaryConstants(Uup,Pup,innerBoundaryUF,innerBoundaryPF,step);
		bigU = thetaVector(iterCounter)*transf2qD(U,x,y,zeroX,zeroY);
		   
		bigUTimeDerivative = -c * domainUtils1Der.YDerivativeZeroBnd(bigU)/h;
	elseif( prmtrs.useZeroBoundary == 3 )

		RD = U(2:end,:);
		RD = flipud(RD);
		bigU = [RD; U];
		
		bigUTimeDerivative = -c * domainUtils1Der.YDerivativeEvenFunZeroBnd(bigU)/h;
        muU = 0;
        muP = 0;
    end
    mu = struct('muU',{muU},'muP',{muP});

    if( prmtrs.useZeroBoundary == 2 )
        [solutionNorms] = CalculateSolutionNorms(U,Uup,P,Pup,UvsUupInfNorm,crrntResidual,...
            residualInfNorm,muU*innerBoundaryUF,muP*innerBoundaryPF,subCounter,thetaVector(iterCounter),step,h);
    else
        solutionNorms = 0;
    end
end

%z1 = YDerivativeEvenFunctions2(P,zeroMatrix,yDerBnd,derivative.second);
%z2 = domainUtils.YDerivativeEvenFunZeroBnd(P);
%mesh(x(zeroX:end),y(zeroY:end),(z2-z1)');
%DrawResidual(x(zeroX:end),y(zeroY:end),al,bt,c,thetaVector(iterCounter),c1,U,derivative.second);
