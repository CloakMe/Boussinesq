clear;clc;
%vpdt2_t10_h025
%vbc_t10_h025

load vpdt2_t10_h025.mat
v = load('vpdt2_t10_h025.dat');
shift = 5; c =2;
[maxvalv, max_err, min_err] = xtrct_prop(x+shift,tt,[1.5 0],v(:,end));

    figure(4)
     for j = 1:size(tt,2) 
        Tr = u_ex(x+5,tt(j),c);% + u_ex(x(:,ff)-5,t(j),-1.5); <-- tochno reshenie
        [maxy maxx] = max(v(:,j));
        %Exact and approximate solution 
        %plot(x, v(:,j),x,Tr,'k');
        %axis([x(1) x(end) -0.03 maxvalv]);
        
        % Oscilations of the solution on a zoomed y-scale [min_err/1 max_err/1]
        plot(x, v(:,j),x,Tr,'k');
        axis([x(1) x(end) min_err/1 max_err/1]);
        
        % The difference b/t exact and approximate solution on a zoomed y-scale [min_err max_err]
        %plot(x, Tr' - v(:,j),x(maxx),0,'ro');
        %axis([x(1) x(end) min_err max_err]);
        F(:,j) = getframe;
        %Error = max(abs(v(:,j)-Tr'))
     end
     movie(F,1)
     Tr = u_ex(x+5,tt(end),c);
     Error = max(abs(Tr' - v(:,end)))