function x=diag3solv(A,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pentsolv.m
%
% Solve a pentadiagonal system Dx=b where D is a strongly nonsingular
% matrix, D is symmetric and each subdiagonal has the following vector form:
% a*ones(1,N-l) == [a a .... a], (l=0,1,2  depending on the diagonal).
% The main diagonal has the following form: [a11 b b .... b a11]
% HERE A = D(2,1:3) !!!!!!!!!!!!!!!!!!!!!!!!
% If D is not a pentadiagonal matrix, results will be wrong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(A,1)*size(A,2)~=3)
    error('A should be a vector, of type [b c b]!');
end
[M,N]=size(b);
if((M>1 && N==1) || (M==1 && N>1))
    if(M>N)
        N=M;
    end
else
    error('b should be a vector, not number or matrix!');
end


x=zeros(N,1);
    
    % Extract bands
    d=A(2);
    f=A(1);
    
    alpha=zeros(N,1);
    gamma=zeros(N-1,1);
    z=zeros(N,1);
    
    % Factor A=LDL'
    alpha(1)=d;
    gamma(1)=f/alpha(1);
    
    alpha(2)=d-f*gamma(1);
    gamma(2)=(f)/alpha(2);
    
    for k=3:N-2
        alpha(k)=d-alpha(k-1)*gamma(k-1)^2;
        gamma(k)=(f)/alpha(k);
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
    
    for k=2:N-1
        x(N-k)=c(N-k)-gamma(N-k)*x(N-k+1);
    end
    
    %for k=N-2:-1:1
    %    x(k)=c(k)-gamma(k)*x(k+1);
    %end
