%clear;clc;
% constants
% st ; end
%-20 ; 15
start_x=-100; end_x = 100;
pw = 0;
h = 0.1/(2)^pw;  tau = 0.05/2^pw;  x = start_x:h:end_x;  t_end=10; 
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
    
    
    %figure(1)
    %plot(x,u_t0,'g',x,dudt_t0,'r',x,u_t02,'co',x,dudt_t02,'k+')
    
ic_utils = IC_2Waves();     
[u00, dudt00] = ic_utils.GetInitialCondition(x, 20);
%t = 10:-0.1:0.0;
t = tt-10;
sol2d = zeros(length(x), length(t));
for p=1:length(t)
    sol2d(:,p) = ic_utils.GetInitialCondition(x, t(p));
end
    
    figure(4)
    mesh(x,t, sol2d')
    xlabel('x')
    ylabel('t')
    colorbar;
    caxis([-.1 .1]);
    view(0,90);
    figure(5)
    mesh(x,t, ((sol2d)-va(:,1:end))')
    colorbar;
    view(0,90);
    xlabel('x')
    ylabel('t')
    
        figure(2)
    plot(x,u00,'g',x,dudt00,'r')
    return;
    k=-1:0.05:1;
    figure(10)
    plot( k, sqrt((k .^ 2 .* (1 - k .^ 2))), 'c')
