function [bigU,bigUTimeDerivative,P,U,bigIC,solutionNorms,...
        theta,c1,c2,iterCounter,zeroX,zeroY,tau,tauVector]=...
        pre_sol_ch_v8(compBox,tau,     h,al,bt1,bt2,c,iterMax,eps,ICSwitch,...
        firstDerivative,secondDerivative,lastP,lastU,lastTheta)

    
x=compBox.x_st:h:compBox.x_end; y=compBox.y_st:h:compBox.y_end;
[zeroX,zeroY]=GetZeroNodes(x,y);

% if(nargin ==15)
%   bigIC = 0;
%   [bigU,bigUTimeDerivative,P,U,theta,c1,c2,iterCounter,tau,solutionNorms,tauVector] =...
%        sol_ch_v8(x,y,lastU,0,iter_max,tau,lastTheta,eps,bt1,bt2,al,c,zeroX,zeroY,lastP);
% 
%    return;
% end
  
    correct = IsExtendedDomainCorrect(compBox);
    bt = bt1/bt2;
    gamma1 =  bt*(1 - c^2 )/(1 - bt*c^2);
    gamma2 =  al*bt/(1 - bt*c^2);
    [X,Y]=Domain(x,y);
    
    
    if(ICSwitch == 1)
       [bigIC,sIC] = pola2cart_v5(x, y, gamma1, gamma2, X, Y);%, ode_name, IC);
    else  
       bigIC = u_ex2d_mat_vc(X,Y,c,bt);
       sIC = bigIC(zeroX:end,zeroY:end);
    end
    th = abs(bigIC(zeroX,zeroY));
    bigIC=bigIC/th;
    sIC = sIC/th;
   
    [bigU,bigUTimeDerivative,P,U,theta,c1,c2,iterCounter,tau,solutionNorms,tauVector] =...
        sol_ch_v8(x,y,sIC    ,bigIC,firstDerivative,secondDerivative,iterMax,tau,th    ,eps,bt1,bt2,al,c,zeroX,zeroY);

if(iterCounter<iterMax && correct == 1 && size(bigUTimeDerivative,1)~=1)
    
    clear('x');clear('y');
    x=compBox.x_st2:h:compBox.x_end2; y=compBox.y_st2:h:compBox.y_end2;
    sx = length(x);
    sy = length(y);
    clear('X');clear('Y');clear('ox1');clear('oy1');
    [X,Y]=Domain(x,y);
    
    clear('bU2');clear('sU2');clear('points');
    c12 = 1-c^2;
    bU2 = c1*theta(iterCounter-1)*(c12*X.^2-Y.^2)./(c12*X.^2+Y.^2).^2;
    points = find_old_grid(compBox.x_st,compBox.x_end,compBox.y_st,compBox.y_end,x,y);
    bU2(points(1):points(2),points(3):points(4)) = bigU*theta(iterCounter-1);
       
    [zeroX,zeroY]=MidPoint(x,y);
       %iter_max = 60000;
       sU2 = bU2(zeroX:end,zeroY:end);
       thet = max(max(abs(bU2)));
       bU2 = bU2/thet;
       sU2 = sU2/thet;
       tau = 1/10000;
       eps = 10^(-8);
       tic
       [bigU,bigUTimeDerivative,P,U,theta,c1,c2,iterCounter,tau,solutionNorms,tauVector] =...
           sol_ch_v8(x,y,sU2,bU2,firstDerivative,secondDerivative,iterMax,tau,thet  ,eps,bt1,bt2,al,c,zeroX,zeroY);
       
end



end
