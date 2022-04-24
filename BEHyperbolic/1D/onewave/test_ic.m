clear;clc;
% constants
% st ; end
%-20 ; 15
start_x=-80; end_x = 120;
pw = 0;
h = 0.05/(2)^pw;  tau = 0.05/2^pw;  x = start_x:h:end_x;  t_end=10; 
beta1=1.5;   beta2=0.5;  alpha=3; beta=beta1/beta2;
%sgm=(4-h^2/tau^2)/12;

c=0;   shift = 5;

estep = max(floor((1/tau)/20),1); %zapazwat se 20 stypki za edinitsa vreme

    %if beta1, beta2 and alpha not specified then beta1 = 1.5; beta2 = 0.5, alpha = 3.
    
    %Tochno reshenie u_ex(x,t,c,alpha,beta1,beta2)
    % t == 0
    u_t0 = u_ex(x+shift,0,c,alpha,beta1,beta2);% + u_ex(x-shift,0,-1.5);
    %Tochno reshenie - proizwodna - dudt_ex(x,t,c,beta1,beta2,alpha)
    % t == 0
    dudt_t0 =dudt_ex(x+shift,0,c,alpha,beta1,beta2);% + dudt_ex(x(l)-shift,0,-1.5);

    %Tochno reshenie u_ex(x,t,c,alpha,beta1,beta2)
    % t == 0
    u_t02 = u_ex(x+shift,35,c,alpha,beta1,beta2);% + u_ex(x-shift,0,-1.5);
    %Tochno reshenie - proizwodna - dudt_ex(x,t,c,beta1,beta2,alpha)
    % t == 0
    dudt_t02 =dudt_ex(x+shift,35,c,alpha,beta1,beta2);% + dudt_ex(x(l)-shift,0,-1.5);
    
    
    figure(1)
    plot(x,u_t0,'g',x,dudt_t0,'r',x,u_t02,'co',x,dudt_t02,'k+')
