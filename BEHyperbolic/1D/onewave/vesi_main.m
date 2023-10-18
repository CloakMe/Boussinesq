clear;clc;
%h = 0.1;  tau = 0.00001;
%''Good" Boussinesq Equation
% with periodic boundary conditions
L1=-100;L2=100;
T0=-15; 
T=20;  
h=0.1; tau=0.01;
alpha = -3;
beta=1;
sigma=0.75;


K1=sqrt(1/8);
K2=-K1;
k=K1;
a1=1;
a2=1;
a12=1;
b=0;
b1=0;
b2=0;

% if abs(K1)>=1||abs(K2)>=1
%     disp('choose |K1|<1 and |K2|<1')
% end
% if K1==-K2
%     disp('choose another K1 and K2')
% end
% 
% if (4*(K1*K1+K1*K2+K2*K2)-3)==0
%     disp('choose another K1 and K2')
% end
%  omega1=-sqrt(K1*K1*(1-K1*K1));
%  omega2=sqrt(K2*K2*(1-K2*K2));

M = (abs(T-T0))/tau;
x = (L1:h:L2)';
N=length(x);

t_end = T - T0;
estep = max(floor((1/tau)/10),1);
e = 1;
tt = zeros(1,t_end*estep);
va = zeros(N,t_end*estep);

%nachalni danni
u0 = vesi_GetInitialCondition(x, T0, k,b, a1, a2, a12)';
u1 = vesi_GetInitialCondition(x, T0+tau, k,b, a1, a2, a12)';
U=u0*0;
figure(1)
plot(x,u0, 'r', x, u1, 'g')

greshka=zeros(floor(M),1);
tic
   
[A,detA]=vesi_create_matrix(N,sigma,h,tau);

 for j = 1:M % sloeve po vremeto
     
   rightside=(2*u1-u0-2*sigma*tau*tau*vesi_deltah(u1,h)+sigma*tau*tau*vesi_deltah(u0,h)+2*sigma*tau*tau*vesi_delta2h(u1,h)-sigma*tau*tau*vesi_delta2h(u0,h)+tau*tau*(vesi_deltah(u1,h)-vesi_delta2h(u1,h)+alpha*vesi_deltah(u1.*u1,h)))';
    U(1:end)=A\rightside;       
      %  U(1:end)=(2*u1-u0+tau*tau*(deltah(u1,h)-delta2h(u1,h)+3*deltah(F,h)))';
        TR = vesi_GetInitialCondition(x, T0+(j+1)*tau, k,b, a1, a2, a12)';
    greshka(j)=max(abs(U-TR));
    u0=u1;
    u1=U;
    if(j==M)
        breakHere = 1
    end
    if(mod(j,estep)==0)
        tt(e)=j*tau;
        va(:,e) = U;
        if(abs(tt(e)-t_end) <= 2*tau)
            TEND = tt(e)
        end
        e=e+1;
    end
end

psi=max(greshka);  
toc

    %plot(x,U,'r',x,TR,'k')
    figure(3)
    mesh(x, tt, va')%(:,1:end-1)
    xlabel('x')
    ylabel('t')
    colorbar;
    view(0,90);

    figure(5)
    plot(x, U, 'r', x, TR, 'k')
          
                      
       