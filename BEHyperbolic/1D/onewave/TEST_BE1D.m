clear;clc;
tic;
% constants

start_x=-60; end_x = 60;
pw = 0;
%h=.2   .0005
h = 0.4;  tau = 0.01;
%h=0.1; tau=0.01;
%h = 0.4;  tau = 0.0005;  
x = start_x:h:end_x;
t_start = 0;
t_interval=10.0;
t_end = 10.0;

beta1=1;   beta2=1;  alpha=-3; beta=beta1/beta2;
sgm = 0.75;
%sgm=(1-h^2/tau^2)/12;

c=3;   shift = 30;
vers = 3;
estep = max(floor((1/tau)/200),1); %zapazwat se 20 stypki za edinitsa vreme
ic_new = true;
if(ic_new)
    ic_utils = IC_2Waves();
    k1 = 1/3;
    k2 = -1/2;
    a1 = 1; 
    a2 = 1; 
    a12 = 1;
    
    %t_start = -20; 
    %b1 = 0;
    %b2 = 0;
    %[u_t0, dudt_t0] = ic_utils.GetInitialCondition2w(x, t_start, k1, k2, b1, b2, a1, a2, a12);
    b1 = -sqrt((k1 .^ 2 .* (1 - k1 ^ 2)))*20;
    b2 = -sqrt((k2 .^ 2 .* (1 - k2 ^ 2)))*20;
    [u_t0, dudt_t0] = ic_utils.GetInitialCondition2w(x, t_start, k1, k2, b1, b2, a1, a2, a12);
    if(nargin == 3)
        
    end
else
    %if beta1, beta2 and alpha not specified then beta1 = 1.5; beta2 = 0.5, alpha = 3.
    
    %Tochno reshenie u_ex(x,t,c,alpha,beta1,beta2)
    % t == 0
     u_t0 = u_ex(x+shift,0,c,alpha,beta1,beta2);
    % u_t0 = u_ex(x+shift,0,c,alpha,beta1,beta2)+u_ex(x-shift,0,-c,alpha,beta1,beta2); 
    %Tochno reshenie - proizwodna - dudt_ex(x,t,c,beta1,beta2,alpha)
    % t == 0
     dudt_t0 =dudt_ex(x+shift,0,c,alpha,beta1,beta2); 
    % dudt_t0 = dudt_ex(x+shift,0,c,alpha,beta1,beta2)+dudt_ex(x-shift,0,-c,alpha,beta1,beta2);
end


    figure(1);plot(x,u_t0,'g',x,dudt_t0,'b');
    title('Initial Condition - u,dudt');
    pause(.1);
    gg=0;
    
%==========================================================================================   
% vz - v-zero layer
% NM --> g_0; O(tau^2 + h^2); 
%  g_1 = alpha*beta*(vz(x_(i-1))^2 + 10*vz(x_i)^2 + vz(x_(i+1))^2 )/12 +  (beta-1)*(vz(x_(i-1) + 10*vz(x_i) + vz(x_(i+1) )/12      
%[v,tt,E,II] = BE1D_v10(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   
% vz - v-zero layer
% NM  O(tau^2 + h^2); numerov along x-axis --> nonlinear term - v^2 only :
%   alpha*beta*(vz(x_(i-1))^2 + 10*vz(x_i)^2 + vz(x_(i+1))^2 )/12      
   %[v,tt,E,II] = BE1D_v11(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   
% IM (g_2)-->Numerov app: g_2 = [ g_0(v_up) + 10*g_0(v_zero) + g_0(v_down) ] / 12  -  if sgm = (1 - h^2/tau^2)/12 ==>  O(tau^4 + h^2) 
  %[v,tt,E,II] = BE1D_vnum1(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   
% IM --> g_1 = [ G(v_up) - G(v_down) ]. / (v_up - v_down)  -  if sgm = (4 + h^2/tau^2)/12 ==>  O(tau^2 + h^2) 
   %[v,tt,va,E,II] = BE1D_v3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   %name = 'v3_O(tau^2 + h^2)_';
   %vers = 1;
   
% IM --> g_3 = 2 * [ G( (1/2)(v_up+v_mid) ) - G( (1/2)(v_mid+v_down) ) ]. / (v_up - v_down) 
% if sgm = (4 + h^2/tau^2)/12 ==>  O(tau^2 + h^4)
   %[vu,tt,va,E,II] = BE1D_v4(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
   
% Taylor  -  O(tau^4 + h^2)
    %[v,dtv,tt,E,II] = BE1D_taylor(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
    %name = 'Taylor_v1_O(tau^4 + h^2)_';
    
% Taylor v2  -   O(tau^4 + h^4)
    %dh = [1 -2 1]/h^2 se zamenq s dh = [-1 16 -30 16 -1]/(12*h^2) i
    %podobrqwame reda na sxodimost w prostranstwenite koordinati
    %[v,dtv,va,tt,II] = BE1D_tv2(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0);
    %name = 'Taylor_v2_O(tau^4 + h^4)_';

% Taylor v3  -   O(tau^4 + h^4)/O(tau^4 + h^2) -> Vasil Vassilev Equation without u_xxtt derivative
    %dh = [1 -2 1]/h^2 se zamenq s dh = [-1 16 -30 16 -1]/(12*h^2) i
    %podobrqwame reda na sxodimost w prostranstwenite koordinati
    [v,dtv,va,tt,II] = BE1D_tv3(start_x,end_x,h,tau,sgm,t_interval,beta1,beta2,alpha,estep,u_t0,dudt_t0,4,4);
    name = 'Taylor_v3_org_O(tau^4 + h^4)_';
    
% Vesi  -   O(tau^2 + h^2)
    %dh = [1 -2 1]/h^2 se zamenq s dh = [-1 16 -30 16 -1]/(12*h^2) i
    %[u_t1, dudt_t1] = ic_utils.GetInitialCondition2w(x, t_start + tau, k1, k2, b1, b2, a1, a2, a12);
    %[v,dtv,va,tt,II, E] = BE1D_vesi(start_x,end_x,h,tau,sgm,t_start + t_interval,beta1,beta2,alpha,estep,u_t0,u_t1, t_start);
    %name = 'vesi_O(tau^2 + h^2)_';
    %vers = 1;
    
    fprintf('elapsed time = %d min\n', toc/60.0);
    save (['SavedWorkspaces\Sol_' name num2str(floor(end_x)) '_tau' num2str(tau * 1000000,'%.07d') '_h0' num2str(h * 100,'%.02d') ]);
%==========================================================================================    
    %figure(2)
    %uex2 =  u_ex(x+shift,tt(end),c,alpha,beta1,beta2);% + u_ex(x(l)-5,t_end,-1.5);
    %plot(x,u_t0,'c',x,v(:,end),'r',x,uex2,'k');
    
    %Error = max(abs(uex2' - v(:,end)))
    
    if(size(II,2)*size(II,1) > 6)
        if(vers == 1)
            Energy = E(end)
            figure(10)
            %plot(t(2:estep:end-1),E,'g',t(2:estep:end-1),E2,'b');
            plot(tt,E-E(1),'c');
            title('The Discrete Energy Functional ');
            xlabel('time t={t}_j_=_1_._.^T'); ylabel('{E}_j_=_1_._.^T');
        end
        Integ = (II(1)+II(end))/2
        
        figure(11)
        plot(tt,II - II(1),'c');
        title('The Discrete Integral Over "u"');
        xlabel('time t={t}_j_=_1_._.^T'); ylabel('{I}_j_=_1_._.^T');
    end

    figure(3)
    mesh(x, tt, va(:,1:length(tt))')
    xlabel('x')
    ylabel('t')
    colorbar;
    view(0,90);
    figure(5)
    plot(x, v, 'g');
    
    %title('End solution');
    if(ic_new)
        hold on;
        u_end = ic_utils.GetInitialCondition2w(x, t_start+t_interval, k1, k2, b1, b2, a1, a2, a12)';
        plot(x(1:10:end), u_end(1:10:end), 'bo');
        hold off;
        set(gca,'FontSize',14);
        legend('numerical solution', 'exact solution');
        fprintf('|u_{ex}-u_h|_L2 = %.10f\n', h*norm(v-u_end,2));
        fprintf('|u_{ex}-u_h|_Inf = %.10f\n', norm(v-u_end,Inf));
    end

    return;
    %vesi
    %|u_{ex}-u_h|_L2 =  0.000131
    %|u_{ex}-u_h|_Inf = 0.000028
    %|u_{ex}-u_h|_L2 =  0.000020
    %|u_{ex}-u_h|_Inf = 0.000003
    %taylor Oh2
    %|u_{ex}-u_h|_L2 =  0.0001390706
    %|u_{ex}-u_h|_Inf = 0.0001499872
    %|u_{ex}-u_h|_L2 =  0.0000240849
    %|u_{ex}-u_h|_Inf = 0.0000367083
    %Oh4
    %|u_{ex}-u_h|_L2 =  0.0000085258
    %|u_{ex}-u_h|_Inf = 0.0000061886
    %|u_{ex}-u_h|_L2 =  0.0000010984
    %|u_{ex}-u_h|_Inf = 0.0000010516
    %movie
    [maxvalv, max_err, min_err] = xtrct_prop(x+shift,tt,u_t0,v(:,end));
   % maxvalv = 3;

    figure(4)
     for j = 1:size(tt,2)
        %Tr = u_ex(x+shift,tt(j),c,alpha,beta1,beta2);% + u_ex(x(:,ff)-5,t(j),-1.5); <-- tochno reshenie
        %[maxy maxx] = max(v(:,j));
        %Approximate solution 
        plot(x, va(:,j));
        %axis([start_x end_x -maxvalv (2*maxvalv)]);
        title('Moving Solution');
        xlabel('x'); ylabel('y');
        % Oscilations of the solution on a zoomed y-scale [min_err/1 max_err/1]
        %plot(x, v(:,j),x,Tr,'k');
        %axis([start_x end_x min_err/1 max_err/1]);
        
        % Oscilations of the solution on a zoomed y-scale [min_err/1 max_err/1]
        %plot(x, v(:,j),x,Tr,'k');
        %axis([x(maxx)-1 x(maxx)+1 1.5-0.1 1.5+0.1]);
        
        % The difference b/t exact and approximate solution on a zoomed y-scale [min_err max_err]
        %plot(x, Tr' - v(:,j),x(maxx),0,'r+',x, Tr*max_err/maxvalv);
        %axis([start_x end_x min_err max_err]);
        
        F(:,j) = getframe;
        %Error = max(abs(v(:,j)-Tr'))
     end
     figure(4)
     movie(F,1,30)
     return;
     %Possibility to save the soluiton:
     %dlmwrite('vta_t10_ht06125.dat', v)
     %dlmwrite('dtvta_t10_h025_t025.dat', dtv)
     %v3d=v(:,end);save v3d.mat;
     %
    sxx = length(x(2601:20:6441));
    x2 = x(1:20:8001);
    j=1;
    figure(13)
    waterfall(x(1:5:end),tt([1 9:10:199]),v(1:5:end,[1 9:10:199])' );
    xlabel('x');   ylabel('time "t"');  zlabel('v');
   
     for j=10:10:200
        plot3(x2,ones(sxx)*tt(j),v(2601:20:6441,j) );
     end
     hold off;
      camva('manual');
      axis('tight')
      view(0,82);
      
      xcon = [ 1e-07,1e-06,1e-05,1e-04,1e-03, 1e-02];
      ycon = [ 3.97, 3.98, 4.16, 3.28, 0.011 0.011];
      yerr1 = [0.0001119103,     0.0001118226,    0.0001109456, 0.0001021745, 0.0001555689, 0.0010459783];
      yerr2 = [0.0000071520,     0.0000070642,    0.000006186, 0.0000105484, 0.000104696, 55];
      
      figure(20)
%       subplot(2,1,1);
%       semilogx(xcon,ycon,'bo')
%       
%         xlabel('\tau');
%         ylabel('Conv speed');
%       subplot(2,1,2);  
      loglog(xcon,yerr1,'d','MarkerFaceColor','r');
      hold on;
      loglog(xcon,yerr2,'kd','MarkerFaceColor','k');
      hold off;
      legend('|u_{0.8} - u^*_{ex}|_\infty', '|u_{0.4} - u^*_{ex}|_\infty');
      text(7e-08, yerr1(1) + 4e-4, '3.97', 'Color', 'green','FontSize',14);
      text(7e-07, yerr1(2) + 4e-4, '3.98', 'Color', 'green','FontSize',14);
      text(7e-06, yerr1(3) + 4e-4, '4.16', 'Color', 'green','FontSize',14);
      text(7e-05, yerr1(4) + 4e-4, '3.28', 'Color', 'green','FontSize',14);
      text(7e-04, yerr1(5) + 4e-4, '0.01', 'Color', 'green','FontSize',14);
      text(3e-03, yerr1(6) + 4e-4, 'none', 'Color', 'green','FontSize',14);
        xlabel('\tau','FontSize',14);
        ylabel('Errors in L_\infty norms','FontSize',14);
        set(gca,'FontSize',14);
        %|u_{ex}-u_h|_L2 = 0.0002339650
%|u_{ex}-u_h|_Inf = 0.0001055689