function [r,dR,R,MAT,dVdr,u]=nat_42(xx,yy,q_end,M,gama1,gama2, X, Y)
    gama2=abs(gama2);
    if(sqrt(xx(end)^2 + yy(end)^2) > q_end)
        r =  0;
        dR = 0;
        R = 0;
        error('r-interval to be evaluated is too large! (r(end) > M');

    end

    sy = length(yy);
    sx = length(xx);
  
    MAT = sqrt(X.^2 + Y.^2);
    vec = MAT(:);
    vec1 = unique(vec);
    vec1 = sort(vec1);
    r = vec1;
    dR = zeros(size(r,1),size(r,2));
    R = dR;
    
    %reshavane na ...
    %q_end=25; %u(x)=c*exp(-x) for x>=q_end
    epsi=10^(-7);
    %M=10000;


    gama=sqrt(gama1);
    h=1/M; %discretization step
    equiDistR=h:h:q_end;
    dVdr=0*equiDistR; %derivative of the solution
    V=0*equiDistR; %the solution

    sizeEquiDistR=size(equiDistR,2);

    %c=1.83;

    a0=gama2*h^3/4;
    h2=h^2;
    b0=1-gama1*h2/2;

    %iteration parameters
    vv=1;
    c_p=0;
    c_n=0;
    cc=1;
    numberOfIter=0;
    %iterations
    while (((vv<-epsi)||(vv>0))&&(numberOfIter<100))

        cas=0;
        %initial conditions (at "infinity")
        dVdr(sizeEquiDistR)=-cc*equiDistR(sizeEquiDistR)*gama*exp(-equiDistR(sizeEquiDistR)*gama);
        V(sizeEquiDistR)=-cc*(1+gama*h/2)*exp(-equiDistR(sizeEquiDistR)*gama);

        for n=sizeEquiDistR-1:-1:1
            %solving for v(n)
            I1=V(n+1);
            a=a0/equiDistR(n);
            b=-b0+gama2*h2*I1;
            c=dVdr(n+1)+h*equiDistR(n)*I1*(gama1+gama2*I1);
            dVdr(n)=2*c/(-b+sqrt(b^2-4*a*c));
            %end solving for v(n)
            if (dVdr(n)>0)
                cas=1;
                break;
            end
            V(n)=V(n+1)+dVdr(n)*h/equiDistR(n);
        end

        if (cas==1)
            c_p=cc;
            if (c_n==0)
                cc=cc/2;
            else
                cc=(c_p+c_n)/2;
            end

        else
            c_n=cc;
            vv=dVdr(1)*M;

            if (c_p==0)
                cc=cc*2;
            else
                cc=(c_p+c_n)/2;
            end
        end

        numberOfIter=numberOfIter+1;

    end

    u=-V;
    cnt = length(r);
    for n=sizeEquiDistR-1:-1:1
        while(equiDistR(n) <= r(cnt)  && r(cnt) < equiDistR(n+1)  )
            dR(cnt) = ( dVdr(n+1) * (r(cnt) - equiDistR(n)) + dVdr(n)* ( equiDistR(n+1) - r(cnt)) )/(equiDistR(n+1) - equiDistR(n));
            R(cnt) = ( u(n+1) * (r(cnt) - equiDistR(n)) + u(n)* ( equiDistR(n+1) - r(cnt)) )/(equiDistR(n+1) - equiDistR(n));
            cnt = cnt - 1;
        end
    end
    dR(1) = dVdr(1);  R(1) = u(1);
        
    N1=sum(equiDistR.*abs(u).^2+abs(dVdr).^2./equiDistR)*h*pi/2;
    N3=sum(equiDistR.*abs(u).^3)*h*pi/2;
    D=N1/2-N3/3;
end
%{
figure(1)
plot(x,v,'y',r,dR,'r')
figure(2)
plot(x,u,'y',r,R,'r')
%}
%{
figure(1)
plot(r,dR,'r')
figure(2)
plot(r,R,'r')
%}

