clear;clc;
% constants

st_x=-5; nd_x = 5;
pw = 0;
h = 0.1;   x = start_x:h:end_x; 
sx = length(x)
b = zeros(length(x),1);
b(1)=0.0005;
b(end)=b(1);
dh=deltah(sx);
dx=deltax(sx);
I = eye(7);
%Idh =  h^2*I+dh;
%sIdh = Idh(2,1:3);
%v = diag3solv(sIdh,b);
%plot(x,v)


v = sech(x);
max(2*(   diag( dh*v')  + diag(dx*v')*dh )-dh*diag(v))
