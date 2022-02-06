clear;clc;
h=0.5; tau = 0.5;
n = 100;
T1 = inv(eye(n)/tau^2 - deltah(n)/(12*h^2));
T2 = inv(h^2*eye(n) - deltah(n));
T3 = deltah(n)*(h^2*eye(n)-deltah(n))/h^2;
n1 = norm(T1,inf)
n2 = norm(T2,inf)
TTT = T2*T1*T3;
DIF = TTT - T1*deltah(n)/h^2;

TTDh = T2*T1*deltah(n)/h^2;

n12 = n1*n2

n3 =norm( T2*T1,inf)
nT = norm(TTT,inf)
nTDh = norm(TTDh,inf)
