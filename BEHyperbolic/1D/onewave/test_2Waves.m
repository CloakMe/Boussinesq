clear;clc;
% constants

start_x=-60; end_x = 60;
pw = 0;
h = 0.1;  tau = 0.01;  x = start_x:h:end_x;  t_end=2.4;
beta1=3;   beta2=1;  alpha=-1; beta=beta1/beta2;
sgm = 1/2;
%sgm=(1-h^2/tau^2)/12;

c=3;   shift = 30;
vers = 3;
estep = max(floor((1/tau)/10),1); %zapazwat se 20 stypki za edinitsa vreme

%if beta1, beta2 and alpha not specified then beta1 = 1.5; beta2 = 0.5, alpha = 3.

%Tochno reshenie u_ex(x,t,c,alpha,beta1,beta2)
% t == 0
% u_t0 = u_ex(x+shift,0,c,alpha,beta1,beta2);
u_t0 = u_ex(x+shift,0,c,alpha,beta1,beta2)+u_ex(x-shift,0,-c,alpha,beta1,beta2); 
%Tochno reshenie - proizwodna - dudt_ex(x,t,c,beta1,beta2,alpha)
% t == 0
% dudt_t0 =dudt_ex(x+shift,0,c,alpha,beta1,beta2); 
dudt_t0 = dudt_ex(x+shift,0,c,alpha,beta1,beta2)+dudt_ex(x-shift,0,-c,alpha,beta1,beta2);
figure(1);plot(x,u_t0,'g',x,dudt_t0,'b');
title('Initial Condition - u,dudt');
pause(.1);

[vUp,tt,v,E,II,vz] = BE1D_v3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,u_t0,dudt_t0); 
vers = 1;
figure(5)
plot(x, vUp, 'r', x, vz, 'b');
title('End solution');  

counter = 0;
while(counter < 300)
    [vUp,tt,va,E,II,vz] = BE1D_v3(start_x,end_x,h,tau,sgm,t_end,beta1,beta2,alpha,estep,vz',dudt_t0,vUp'); 
    if(size(II,2)*size(II,1) > 6)
        if(vers == 1)
            Energy = E(end)
            figure(10)
            %plot(t(2:estep:end-1),E,'g',t(2:estep:end-1),E2,'b');
            plot(tt,E,'c');
            title('The Discrete Energy Functional ');
            xlabel('time t={t}_j_=_1_._.^T'); ylabel('{E}_j_=_1_._.^T');
        end
        Integ = (II(1)+II(end))/2

        figure(11)
        plot(tt,II,'c');
        title('The Discrete Integral Over "u"');
        xlabel('time t={t}_j_=_1_._.^T'); ylabel('{I}_j_=_1_._.^T');
    end

    %movie
    [maxvalv, max_err, min_err] = xtrct_prop(x+shift,tt,u_t0,vUp(:,end));
    % maxvalv = 3;
    figure(5)
    plot(x, vUp, 'r', x, vz, 'b');
    title('End solution');
    counter = counter+1;
    temp = vz;
    vz = vUp;
    vUp = temp;
end

figure(2)
plot(x, u_t0, 'r', x, vUp, 'b');