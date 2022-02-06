function x=pentsolv(a11,A,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pentsolv.m
%
% Solve a pentadiagonal system Dx=b where D is a strongly nonsingular
% matrix, D is symmetric and each subdiagonal has the following vector form:
% a*ones(1,N-l) == [a a .... a], (l=0,1,2  depending on the diagonal).
% The main diagonal has the following form: [a11 b b .... b a11]
% HERE A = D(3,1:5) !!!!!!!!!!!!!!!!!!!!!!!!
% If D is not a pentadiagonal matrix, results will be wrong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(A,1)*size(A,2)~=5)
    error('A should be a vector, of type [a b c b a]!');
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
    d=A(3);
    f=A(2);
    e=A(1);
        
    alpha=zeros(N,1);
    gamma=zeros(N-1,1);
    delta=zeros(N-2,1);
    z=zeros(N,1);
    
    % Factor A=LDL'
    alpha(1)=a11;
    gamma(1)=f/alpha(1);
    delta(1)=e/alpha(1);
    
    alpha(2)=d-f*gamma(1);
    gamma(2)=(f-e*gamma(1))/alpha(2);
    delta(2)=e/alpha(2);
    
    for k=3:N-2
        alpha(k)=d-e*delta(k-2)-alpha(k-1)*gamma(k-1)^2;
        gamma(k)=(f-e*gamma(k-1))/alpha(k);
        delta(k)=e/alpha(k);
    end
    
    alpha(N-1)=d-e*delta(N-3)-alpha(N-2)*gamma(N-2)^2;
    gamma(N-1)=(f-e*gamma(N-2))/alpha(N-1);
    alpha(N)=a11-e*delta(N-2)-alpha(N-1)*gamma(N-1)^2;
    
    % Update Lx=b, Dc=z
    
    z(1)=b(1);
    z(2)=b(2)-gamma(1)*z(1);
    
    for k=3:N
        z(k)=b(k)-gamma(k-1)*z(k-1)-delta(k-2)*z(k-2);
    end
    
    c=z./alpha;
        
    % Backsubstitution L'x=c
    x(N)=c(N);
    x(N-1)=c(N-1)-gamma(N-1)*x(N);
    
    for k=2:N-1
        x(N-k)=c(N-k)-gamma(N-k)*x(N-k+1)-delta(N-k)*x(N-k+2);
    end
    