function [bigU,bigUTimeDerivative,P,U,bigIC,solutionNorms,theta,c1,c2,tauVector,angl]=...
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
       bigIC = u_ex2d_mat_vc(X,Y,c/sqrt(bt),bt);
       IC = bigIC(zeroX:end,zeroY:end);
    end
    th = abs(IC(1,1));
    IC = IC/th;
   
    [bigU,bigUTimeDerivative,P,U,theta,c1,c2,solutionNorms,tauVector,angl] =...
        sol_ch_v8(IC,x,y,prmtrs,bt1,bt2,al,c,th,zeroX,zeroY,derivative);
        
end
