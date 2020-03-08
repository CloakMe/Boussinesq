function [x,t,u0] = BPE_symplectic_periodic_condition( L, T, h, tau)
    %Boussinesq Paradigm Equation
    % with periodic boundary conditions
    alfa=3; beta1=1.5; beta2=0.5;
    c=2;
    % mreza
    M = T/tau;
    x = (-L:h:L)';
    N=length(x);
    % nachalni danni
    u0 =exact_sol_p2(x,0,c,beta1,beta2,alfa);
    v0=-c*exact_sol_p2(x,0,c,beta1,beta2,alfa);
    t=zeros(M,1);
    
    % vytreshnite etapi na metoda                    
    k1=0*u0;k2=0*u0; k3=0*u0;                       
    l1=0*v0;l2=0*v0; l3=0*v0;                       

    greshka=zeros(M,1);

    a=-beta1/(4*h*h);   C=(1+2*beta1/(4*h*h));
    A=create_matrix(N,a,C);

    % koeficienti na Partitioned Runge-Kutta method
    A11=0;A12=0;A13=0;
    A21=5/24; A22=1/3; A23=-1/24;
    A31=1/6; A32=2/3; A33=1/6;
    B1=1/6; B2=2/3; B3=1/6;

    a11=1/6; a12=-1/6; a13=0;
    a21=1/6; a22=1/3; a23=0;
    a31=1/6; a32=5/6; a33=0;

    b1=1/6; b2=2/3; b3=1/6;
    maxiter=1;
    for k = 1:M % sloeve po vremeto
        
        k1=xhat(v0,h);	k2=xhat(v0,h); k3=xhat(v0,h);
        vec1=u0+tau*(a11*k1+a12*k2+a13*k3);
        right1=B(vec1,alfa,beta2,h);
        l1=A\right1;

        vec2=u0+tau*(a21*k1+a22*k2+a23*k3);
        right2=B(vec2,alfa,beta2,h);
        l2=A\right2;

        vec3=u0+tau*(a31*k1+a32*k2+a33*k3);
        right3=B(vec3,alfa,beta2,h);
        l3=A\right3;

        %nuleva iteraciq  
        % iteracionen metod
        sub = 1;iter=1;
        while sub > 10^(-14)

            RIGHT1=v0+tau*(A11*l1+A12*l2+A13*l3);
            RIGHT2=v0+tau*(A21*l1+A22*l2+A23*l3);
            RIGHT3=v0+tau*(A31*l1+A32*l2+A33*l3);
            kk1=xhat(RIGHT1,h);	kk2=xhat(RIGHT2,h); kk3=xhat(RIGHT3,h);
            vec11=u0+tau*(a11*kk1+a12*kk2+a13*kk3);
            right11=B(vec11,alfa,beta2,h);
            ll1=A\right11;

            vec22=u0+tau*(a21*kk1+a22*kk2+a23*kk3);
            right22=B(vec22,alfa,beta2,h);
            ll2=A\right22;

            vec33=u0+tau*(a31*kk1+a32*kk2+a33*kk3);
            right33=B(vec33,alfa,beta2,h);
            ll3=A\right33;

            RAZ(1)=max(abs(l1-ll1));RAZ(2)=max(abs(l2-ll2)); RAZ(3)=max(abs(l3-ll3));
            RAZ(4)=max(abs(k1-kk1));RAZ(5)=max(abs(k2-kk2)); RAZ(6)=max(abs(k3-kk3));
            sub=max(RAZ);
            iter=iter+1;
            k1=kk1;k2=kk2;k3=kk3;
            l1=ll1;l2=ll2;l3=ll3;

            if iter>100
                break
            end
        end

        if iter > maxiter
            maxiter=iter;
        end

        U=u0+tau*(b1*k1+b2*k2+b3*k3);
        V=v0+tau*(B1*l1+B2*l2+B3*l3);
        u0=U;
        v0=V;
        TR= exact_sol_p2(x,(k)*tau,c,beta1,beta2,alfa);
        greshka(k)=max(abs(u0-TR));
        t(k+1) = tau*k;
    end
    psi=max(greshka);
end