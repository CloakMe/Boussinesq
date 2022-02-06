clear;clc;
% constants
% st ; end
%-20 ; 15
startx=-20;
tau = 0.05;  h = 0.05;  sgm=1/2;  x = startx:h:25;  t_end=6;
beta1=1.5;   beta2=0.5;  alpha=3;
beta=beta1/beta2;
t = 0:tau:t_end;

sx = size(x,2);

  
figure(1)
     for j = 1:size(t,2) 
        Tr = u_ex(x+5,t(j),2) + u_ex(x-5,t(j),-1.5);
        plot(x,Tr,'k');
        axis([startx x(end,end) (-0.3) 2.5]);
        F(:,j) = getframe;
     end
     movie(F,1)

    