function x=tridiagsolv(A,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pentsolve.m
%
% Solve a pentadiagonal system Ax=b where A is a strongly nonsingular matrix
% 
% If A is not a pentadiagonal matrix, results will be wrong
%
% Reference: G. Engeln-Muellges, F. Uhlig, "Numerical Algorithms with C"
%               Chapter 4. Springer-Verlag Berlin (1996)
%
% Written by Greg von Winckel 3/15/04
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(size(A,1)*size(A,2)~=3)
    error('A should be a vector, of type [b c b]!');
    return;
end
[M,N]=size(b);
if((M>1 && N==1) || (M==1 && N>1))
    if(M>N)
        N=M;
    else
        M=N;
    end
else
    error('b should be a vector, not number or matrix!');
    return;
end

x=zeros(N,1);
    
   % Extract bands
    d=A(2);
    f=A(1);
           
    alpha=zeros(N,1);
    gamma=zeros(N-1,1);
    delta=zeros(N-2,1);
    c=zeros(N,1);
    z=zeros(N,1);
    
    % Factor A=LDL'
    alpha(1)=d;
    gamma(1)=f/alpha(1);
        
    alpha(2)=d-f*gamma(1);
    gamma(2)=f;
      
    for k=3:N-2
        alpha(k)=d-alpha(k-1)*gamma(k-1)^2;
        gamma(k)=f;
    end
    
    alpha(N-1)=d-alpha(N-2)*gamma(N-2)^2;
    gamma(N-1)=(f)/alpha(N-1);
    alpha(N)=d-alpha(N-1)*gamma(N-1)^2;
    
    % Update Lx=b, Dc=z
    
    z(1)=b(1);
    z(2)=b(2)-gamma(1)*z(1);
    
    for k=3:N
        z(k)=b(k)-gamma(k-1)*z(k-1);
    end
    
    c=z./alpha;
        
    % Backsubstitution L'x=c
    x(N)=c(N);
    x(N-1)=c(N-1)-gamma(N-1)*x(N);
    
    for k=N-2:-1:1
        x(k)=c(k)-gamma(k)*x(k+1);
    end
    
     
    