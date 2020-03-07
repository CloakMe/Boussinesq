function A=create_matrix(N,a,c)
    % vryshta matrica ot vida
    % c  0 a 0......0 a 0
    % 0  c 0 a......0 0 a
    %....a 0 c 0 a...........
    %a  0  0 0....a 0 c 0
    %0  a  0 0... 0 a 0 c
    Z=zeros(N,N);
    Z(1,1)=c; Z(1,3)=a; Z(1,N-1)=a;
    Z(2,2)=c; Z(2,4)=a;Z(2,N)=a;
    
    for i=3:N-2
        Z(i,i-2)=a; Z(i,i)=c; Z(i,i+2)=a;
    end
    
    Z(N-1,1)=a; Z(N-1,N-3)=a; Z(N-1,N-1)=c;
    Z(N,2)=a;                            Z(N,N-2)=a; Z(N,N)=c;
    A=sparse(Z);
end