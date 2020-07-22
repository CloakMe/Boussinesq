function [bigU,bigUTimeDerivative,Pup,Uup,thetaVector,mu,solutionNorms,tauVector,angl,sw_div]=...
    sol_ch_v8(U,x,y,prmtrs,bt1,bt2,al,c,theta,derivative,P)

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
    [approxBoundaryF2, approxBoundaryF1] = GetApproximationForBoundary(x(zeroX:end),y(zeroY:end),h,c);

    outerTopBoundaryF1=approxBoundaryF1(1:end-4,end-3:end); 
    outerRigthBoundaryF1=approxBoundaryF1(end-3:end,1:end-4);
    outerTopBoundaryF2=approxBoundaryF2(1:end-4,end-3:end); 
    outerRigthBoundaryF2=approxBoundaryF2(end-3:end,1:end-4);

    xLength=length(x(zeroX:end));  yLength=length(y(zeroY:end));

    [yLen,xLen,yHeiA,xHeiA] = GetInnerNodesForComparingBoundaryFunctions(xLength,yLength,step);

    innerBoundaryUF1=[approxBoundaryF1(xHeiA,yLen)  approxBoundaryF1(xLen,yHeiA)'];
    innerBoundaryPF1 = bt*(1-c^2)*innerBoundaryUF1;
    innerBoundaryUF2=[approxBoundaryF2(xHeiA,yLen)  approxBoundaryF2(xLen,yHeiA)'];
    innerBoundaryPF2 = bt*(1-c^2)*innerBoundaryUF2;

    bigZeroMatrix = zeros(sx,sy);
    zeroMatrix = zeros(size(U));
    if(prmtrs.useZeroBoundary == 0)
        [muU1,muU2,muP1,muP2] = FindBoundaryConstantsExtended(U,...
            0*U,...
            innerBoundaryUF1,...
            innerBoundaryUF2,...
            innerBoundaryPF1,...
            innerBoundaryPF2,...
            step);
    else
        muU1 = 0;
        muU2 = 0;
    end
    deltaU = DeltaEvenFunctions(...
        U,...
        zeroMatrix,...
        muU1*outerRigthBoundaryF1 + muU2*outerRigthBoundaryF2,...
        muU1*outerTopBoundaryF1 + muU2*outerTopBoundaryF2,...
        derivative.second);
    
    if(sw == 0)
        P = bt*(1-c^2)*U - (1-bt*c^2)*deltaU + theta*al*bt*U.^2;
    end
     
     iterCounter=1;
     subCounter=1;
     Usquare = U .^2;
     thetaVector(iterCounter) = (P(1,1) - bt*(1-c^2)*U(1,1) + (1-bt*c^2)*deltaU(1,1) )/(al*bt*Usquare(1,1));
     residualInfNorm(1) =  GetResidualInfNorm(...
         al,...
         bt,...
         c,...
         thetaVector,...
         iterCounter,...
         U,...
         Usquare,...
         deltaU,...
         zeroMatrix,...
         muU1*outerRigthBoundaryF1 + muU2*outerRigthBoundaryF2,...
         muU1*outerTopBoundaryF2 + muU2*outerTopBoundaryF2,...
         derivative.second);
     
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
        if(prmtrs.useZeroBoundary == 0)
            %[muU1,c2] = FindBoundaryConstants(U,P,innerBoundaryUF1,innerBoundaryPF1,step);
            [muU1,muU2,muP1,muP2] = FindBoundaryConstantsExtended(U,...
                P,...
                innerBoundaryUF1,...
                innerBoundaryUF2,...
                innerBoundaryPF1,...
                innerBoundaryPF2,...
                step);
        else
            muU1 = 0; muU2 = 0; muP1 = 0; muP2 = 0;
        end
        outerRigthBoundaryF = muU1*outerRigthBoundaryF1 + muU2*outerRigthBoundaryF2;
        outerTopBoundaryF = muU1*outerTopBoundaryF1 + muU2*outerTopBoundaryF2;
        deltaU = DeltaEvenFunctions(U, zeroMatrix, outerRigthBoundaryF, outerTopBoundaryF, derivative.second); %<--- TIME CONSummable
        Usquare = U .^2;
        thetaVector(iterCounter) = (P(1,1) - bt*(1-c^2)*U(1,1) + (1-bt*c^2)*deltaU(1,1) )/(al*bt*Usquare(1,1));
         % TIME CONSumable       |
         %                      \|/
         %                       '
        yDerBnd = bt*(1-c^2)*(muP1 * outerTopBoundaryF1 + muP2 * outerTopBoundaryF2);
        xDerBnd = bt * c^2 * outerRigthBoundaryF + bt * (1-c^2) * (muP1*outerRigthBoundaryF1 + muP2*outerRigthBoundaryF2);
        Pup = P + (tau)*(YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second) +...
            XDerivativeEvenFunctions(P-bt*c^2*(deltaU - U),zeroMatrix,xDerBnd,derivative.second));  
        
        Uup = U + tau*( Pup - (al*bt*thetaVector(iterCounter) )*Usquare +...
            ((1 - bt*c^2))* deltaU - (bt*(1 - c^2))*U );
        
        UvsUupInfNorm(iterCounter) = max(max(abs(U-Uup)));
        
        if(mod(iterCounter,10) ==0)        
           
           crrntResidual = Get2DResidual(al,bt,c,thetaVector,iterCounter,U,Usquare,deltaU,...
               zeroMatrix,outerRigthBoundaryF,outerTopBoundaryF,derivative.second);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=abs(thetaVector(iterCounter))*max(max(abs(crrntResidual(1:end-8,1:end-8))));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           [flag, Px, Py] = StopCriteria(x, y, zeroX, zeroY, U, ax, ay, minResidual, eps);
           
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               %fprintf('|ax - axNew| = %.15e \n|ay - ayNew| = %.15e\n eps = %.15e\n', abs( ax-Px(1) ), abs( ay-Py(1) ), eps );
                
               if(prmtrs.plotResidual)
                   PlotResidual(x(zeroX:end),y(zeroY:end),crrntResidual*thetaVector(iterCounter));
               end
               if(prmtrs.plotBoundary)
                   PlotBoundary(x,y,zeroX, U, outerTopBoundaryF);
               end
               if(prmtrs.plotAssympt)
                   PlotAssymptVsSolu( x(zeroX:end), y(zeroY:end), h, 1, 1, U, muU2*thetaVector(iterCounter), c);
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
    [muU1, muU2, muP1, muP2]=FindBoundaryConstantsExtended(Uup,Pup,innerBoundaryUF1,innerBoundaryUF2,innerBoundaryPF1,innerBoundaryPF2,step);
    mu = struct('muU1',{muU1},'muU2',{muU2},'muP1',{muP1},'muP2',{muP2});

    [bigU,bigUTimeDerivative] = GetBigSolution(x,y,c,bt,thetaVector(iterCounter),...
        Uup,bigZeroMatrix,derivative.first,outerTopBoundaryF);
    innerBoundaryUF = muU1 * innerBoundaryUF1 + muU2 * innerBoundaryUF2;
	innerBoundaryPF = muP1 * innerBoundaryPF1 + muU2 * innerBoundaryPF2;
    [solutionNorms] = CalculateSolutionNorms(U,Uup,P,Pup,UvsUupInfNorm,crrntResidual,...
        residualInfNorm,innerBoundaryUF,innerBoundaryPF,subCounter,thetaVector(iterCounter),step,h);
end

%z1 = YDerivativeEvenFunctions2(P,zeroMatrix,yDerBnd,derivative.second);
%z2 = YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second);
%mesh(x(zeroX:end),y(zeroY:end),(z2-z1)');
%DrawResidual(x(zeroX:end),y(zeroY:end),al,bt,c,thetaVector(iterCounter),c1,U,derivative.second);
