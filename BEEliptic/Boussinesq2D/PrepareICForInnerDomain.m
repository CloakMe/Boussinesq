function [bigU,bigUTimeDerivative,P,U,bigIC,solutionNorms,theta,mu,tauVector,angl,sw_div]=...
        PrepareICForInnerDomain(compBox,prmtrs,al,bt1,bt2,c,derivative)%,lastP,lastU,lastTheta)
        

    x=compBox.x_st:prmtrs.h:compBox.x_end; y=compBox.y_st:prmtrs.h:compBox.y_end;
    [zeroX,zeroY]=GetZeroNodes(x,y);
  
    bt = bt1/bt2;
    gamma1 =  (bt - c^2 )/(1 - c^2);
    gamma2 =  al*bt/(1 - c^2);
    [X,Y]=Domain(x,y);
    
    if(prmtrs.ICSwitch == 1)
       [bigIC,IC] = pola2cart_v5(x, y, gamma1, gamma2);%, ode_name, IC);
    else
       bigIC = u_ex2d_mat_vc(X,Y,c,bt);
       IC = bigIC(zeroX:end,zeroY:end);
    end
    th = abs(IC(1,1));
    IC = IC/th;
    
    if(prmtrs.useZeroBoundary == 2)
        [bigU,bigUTimeDerivative,P,U,theta,mu,solutionNorms,tauVector,angl,sw_div] =...
            sol_ch_v8ZeroBnd(IC,x,y,prmtrs,bt1,bt2,al,c,th,derivative);
    elseif(prmtrs.useZeroBoundary == 3)
        
		IC_half = bigIC(zeroX:end,:)/th;
        figure(2); mesh(x(zeroX:end),y,IC_half');
        
        [bigU,bigUTimeDerivative,P,U,theta,mu,solutionNorms,tauVector,angl,sw_div] =...
            sol_ch_v8ZeroBnd(IC_half,x,y,prmtrs,bt1,bt2,al,c,th,derivative);
    else
        [bigU,bigUTimeDerivative,P,U,theta,mu,solutionNorms,tauVector,angl,sw_div] =...
            sol_ch_v8(IC,x,y,prmtrs,bt1,bt2,al,c,th,derivative);
    end    
        
end
