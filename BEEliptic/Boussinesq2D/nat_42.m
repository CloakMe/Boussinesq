function [r,dR,R,MAT,v,u]=nat_42(xx,yy,q_end,M,gama1,gama2, X, Y)
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
   cnt = length(r);
   
   %reshavane na ...
   %q_end=25; %u(x)=c*exp(-x) for x>=q_end
   epsi=10^(-7);
   %M=10000;


    gama=sqrt(gama1);
    h=1/M; %discretization step
    x=h:h:q_end;
    v=0*x; %derivative of the solution
    I=0*x; %the solution

    nn=size(x,2);

    %c=1.83;

    a0=gama2*h^3/4;
    h2=h^2;
    b0=1-gama1*h2/2;

    %iteration parameters
    vv=1;
    c_p=0;
    c_n=0;
    cc=1;
    nnn=0;
%iterations
while (((vv<-epsi)||(vv>0))&&(nnn<100))

    cas=0;
    %initial conditions (at "infinity")
    v(nn)=-cc*x(nn)*gama*exp(-x(nn)*gama);
    I(nn)=-cc*(1+gama*h/2)*exp(-x(nn)*gama);

    for n=nn-1:-1:1
    %solving for v(n)
    I1=I(n+1);
    a=a0/x(n);
    b=-b0+gama2*h2*I1;
    c=v(n+1)+h*x(n)*I1*(gama1+gama2*I1);
    v(n)=2*c/(-b+sqrt(b^2-4*a*c));
    %end solving for v(n)
    if (v(n)>0)
        cas=1;
    break;
    end
        I(n)=I(n+1)+v(n)*h/x(n);
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
            vv=v(1)*M;

            if (c_p==0)
                cc=cc*2;
            else
            cc=(c_p+c_n)/2;
            end
        end

        nnn=nnn+1;

 end

   u=-I;
        for n=nn-1:-1:1
            while(x(n) <= r(cnt)  && r(cnt) < x(n+1)  )
                dR(cnt) = ( v(n+1) * (r(cnt) - x(n)) + v(n)* ( x(n+1) - r(cnt)) )/(x(n+1) - x(n));
                R(cnt) = ( u(n+1) * (r(cnt) - x(n)) + u(n)* ( x(n+1) - r(cnt)) )/(x(n+1) - x(n));
                cnt = cnt - 1;
            end
        end
        dR(1) = v(1);  R(1) = u(1);
        
  N1=sum(x.*abs(u).^2+abs(v).^2./x)*h*pi/2;
  N3=sum(x.*abs(u).^3)*h*pi/2;
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

