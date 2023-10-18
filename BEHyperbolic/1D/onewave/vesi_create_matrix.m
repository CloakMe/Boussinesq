function [A,detA]=vesi_create_matrix(N,sigma,h,tau)
% vryshta matrica ot vida
% c      -4a-b    a   0......   0 
% -4a-b   c     -4a-b 0......   0
% a     -4a-b    c    -4a-b 0...0
% a       0      0      0   ....a -4a-b   c  -4a-b
%0  0  0 0               ...    0    a  -4a-b   c
a=sigma*tau*tau/(h*h*h*h);
b=sigma*tau*tau/(h*h);
c=1+2*b+6*a;

Z=zeros(N,N);
Z(1,1)=c; Z(1,2)=-4*a-b; Z(1,3)=a;
Z(2,1)=-4*a-b; Z(2,2)=c;     Z(2,3)=-4*a-b;   Z(2,4)=a;
for i=3:N-2
    Z(i,i-2)=a; Z(i,i-1)=-4*a-b; Z(i,i)=c; Z(i,i+1)=-4*a-b; Z(i,i+2)=a;
end
Z(N-1,N-3)=a; Z(N-1,N-2)=-4*a-b;    Z(N-1,N-1)=c;  Z(N-1,N)=-4*a-b;
                Z(N,N-2)=a;  Z(N,N-1)=-4*a-b; Z(N,N)=c;
                detA=det(Z);
A=sparse(Z);

end