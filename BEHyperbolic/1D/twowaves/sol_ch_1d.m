function [Fup,Uup,residualInfNorm,thetaVector,tauVector]=...
    sol_ch_1d(U,x,prmtrs,bt1,bt2,al,c,derivative)

% constants
    tau = prmtrs.tau;
    h=prmtrs.h;
    eps = prmtrs.eps;
    iterMax = prmtrs.iterMax;
    bt = bt1/bt2;
    
    UvsUupInfNorm = ones(1,iterMax);
    thetaVector = zeros(1,iterMax);
    tauVector = zeros(1,iterMax); 
    tauVector(1) = tau;
    residualInfNorm = ones(1,iterMax/10);

    zeroMatrix = zeros(size(U));

     iterCounter=1;
     subCounter=1;
        
     Uxx = XDerivativeEvenFunctions(U, zeroMatrix,0,derivative.second);
     F = bt*(1-c^2)*U - (1-bt*c^2)*Uxx; 
     nonlinTerm = XDerivativeEvenFunctions(U,zeroMatrix,0,derivative.first).^2 + 2*U.*Uxx;
     thetaVector(iterCounter) = (F(1))/(4*al*bt*U(1)*Uxx(1));
     residualInfNorm(1) = abs(thetaVector(iterCounter))*max(XDerivativeEvenFunctions(F,zeroMatrix,0,derivative.second)-2*al*bt*nonlinTerm );
     
     Fup = F;
     Uup = U;
     FD = derivative.second;
     M = FD(1);
     mid = ( length( derivative.second ) + 1 ) / 2;% ( this.order/2 + 1 );
     if( mid == 2 )
         M = FD(1,1);
     elseif( mid == 3 )   
         M = [ FD(1,1:2); 0, FD(1,1) ];
     elseif( mid == 4 )   
         M = [ FD(1,1:3); 0, FD(1,1:2); 0, 0, FD(1,1) ];
     end
      
     minResidual = min(1000,residualInfNorm(subCounter));

     fig9=figure(9);
     set(gcf,'currentcharacter','T');
     sw_div = 0; stopFlag = 0; manualStop=0;
    while( 1 > 0 ) %~( residualInfNorm(subCounter) < eps )
        F = Fup;
        U = Uup;
        iterCounter=iterCounter+1;

        Uxx = XDerivativeEvenFunctions(U, zeroMatrix,0,derivative.second);
        
        nonlinTerm = XDerivativeEvenFunctions(U,zeroMatrix,0,derivative.first).^2 + 2*U.*Uxx;
        
         % TIME CONSumable       |
         %                      \|/
         %                       '
        thetaVector(iterCounter) = (F(1))/(4*al*bt*U(1)*Uxx(1));
        
        Fup = F + tau*(2*al*bt*thetaVector(iterCounter) * nonlinTerm - XDerivativeEvenFunctions(F,zeroMatrix,0,derivative.second) );  
        
        Uup = U + tau*( Fup - (bt*(1 - c^2))*U - (bt*c^2-1)*Uxx  );        
        
        UvsUupInfNorm(iterCounter) = max(max(abs(U-Uup)));
        
        if(mod(iterCounter,10) ==0)
           
           crrntResidual = thetaVector(iterCounter)*(XDerivativeEvenFunctions(F,zeroMatrix,0,derivative.second)-2*al*bt*nonlinTerm);
           
           subCounter=subCounter+1;
           [residualInfNorm(subCounter)]=max(crrntResidual(1:end-8));
           minResidual = min(minResidual,residualInfNorm(subCounter));

           if(minResidual < eps)
               stopFlag = 1;
           end
           
           if(mod(iterCounter,500) ==0)
               fprintf('%d \n',iterCounter);
               fprintf('||R||_Inf = %.4e \n', residualInfNorm(subCounter));
               fprintf('tau = %.4e \n', tau);
               
               manualStop = Stop(fig9,manualStop);
           end
           
           if( residualInfNorm(subCounter) > 40000 || ...
               isnan(residualInfNorm(subCounter)) || ...
               isinf(residualInfNorm(subCounter)))

               warning(' Algo diverges!');
               sw_div = 1;
           end
           
        end
        
        tauVector(iterCounter) = tau;
               
        if(sw_div ==1 || manualStop ==2 || iterCounter ==iterMax || stopFlag ==1)
            break;
        end
    end
    
   % UvsUupInfNorm = thetaVector(iterCounter)*UvsUupInfNorm(1:iterCounter);
    tauVector = tauVector(1:iterCounter);
    thetaVector = thetaVector(1:iterCounter);
    residualInfNorm = residualInfNorm(1:subCounter);
%     [bigU,bigUTimeDerivative] = GetBigSolution(x,y,c,bt,thetaVector(iterCounter),...
%         Uup,bigZeroMatrix,derivative.first,outerTopBoundaryF*muU);
%     
%     [solutionNorms] = CalculateSolutionNorms(U,Uup,F,Fup,UvsUupInfNorm,crrntResidual,...
%         residualInfNorm,muU*innerBoundaryUF,muP*innerBoundaryPF,subCounter,thetaVector(iterCounter),step,h);
end

%z1 = YDerivativeEvenFunctions2(P,zeroMatrix,yDerBnd,derivative.second);
%z2 = YDerivativeEvenFunctions(P,zeroMatrix,yDerBnd,derivative.second);
%mesh(x(zeroX:end),y(zeroY:end),(z2-z1)');
%DrawResidual(x(zeroX:end),y(zeroY:end),al,bt,c,thetaVector(iterCounter),c1,U,derivative.second);
%mesh(x(zeroX:end),y(zeroY:end),PTilda');
%mesh(x(zeroX:end),y(zeroY:end),Pup');
%mesh(x(zeroX:end),y(zeroY:end),Pup1');
