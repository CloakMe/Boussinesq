function [tau, tauMax, tauIncreasedIteration, tauDecreasedIteration] =...
    DefineCurrentTau(subCounter, iterCounter, iterMax, tau,  tauMax, tauIncreasedIteration,...
    tauDecreasedIteration, rsdlInf, smallUdiffInf, minResidual, angl, stopSwitch, crrntResidual)

    if( iterCounter >  tauDecreasedIteration && iterCounter >  500)
        if( (iterCounter > 75 &&  rsdlInf(subCounter)> 50*minResidual) )
            
            tauDecreasedIteration=iterCounter+10;
            tauMax = 0.98*tau;  
            tau = 0.72*tau; 
        end
        vecX = abs( diff( sign( crrntResidual(floor((1+end)/2):end, 1) ) ) );
        vecY = abs( diff( sign( crrntResidual(1, floor((1+end)/2):end) ) ) );
        vec = histc(vecX, unique(vecX)) + histc(vecX, unique(vecY));
        flag = (abs(crrntResidual(1,1) - crrntResidual(1,3))*3 < abs(crrntResidual(1,1) - crrntResidual(1,2)) &&...
            abs(crrntResidual(1,2) - crrntResidual(1,4))*3 < abs(crrntResidual(1,2) - crrntResidual(1,3)) );
        if(iterCounter > 75 && flag == 1)
           %if(vec(2)/sum(vec) > 0.5) 
               tauDecreasedIteration=iterCounter+10;
               %tauMax = 0.98*tau; 
               tau = 0.95*tau; 
           %end
        end
%         if( (iterCounter > 75 &&  abs(angl(subCounter)) > 180 + tau ) )
%             tauDecreasedIteration=iterCounter+10;
%             tau = 0.95*tau; 
%         end
        %exponent = 10^(-ceil(log10(rsdlInf(iterCounter-50)))+1);
        %minResidualCurrent = rsdlInf(iterCounter)*exponent;
        %wtf = rsdlInf(iterCounter-50)*exponent- minResidualCurrent;
        
        if( smallUdiffInf(iterCounter-10)>smallUdiffInf(iterCounter) &&...
            angl(subCounter) < 180 &&  tau < tauMax &&...  % && 0<wtf && wtf<(0.3+minResidualCurrent/30)
          iterCounter>tauIncreasedIteration+10 && iterCounter<iterMax - 1400 &&...
          stopSwitch==0 && rsdlInf(subCounter)<= minResidual*2)
      
            tauIncreasedIteration = iterCounter;
            tau = 1.015*tau; 
        end
    end

end

