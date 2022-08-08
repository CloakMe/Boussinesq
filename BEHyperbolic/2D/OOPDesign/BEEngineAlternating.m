classdef (ConstructOnLoad) BEEngineAlternating < BEEngine
   % Class help goes here
  
  methods
    
    function this = BEEngineAlternating( dscrtParams, eqParams, ic, flag )
        if(nargin == 3)
            flag = 0;
        end
        this = this@BEEngine( dscrtParams, eqParams, ic, flag );
        
        sy = length(this.y);
        sx = length(this.x);
        longerSideNodes = (1+(max(sx,sy)-min(sx,sy))/2):(max(sx,sy)-(max(sx,sy)-min(sx,sy))/2);
        if(length(this.x) > length(this.y))
            xnodes = longerSideNodes;
            ynodes = 1:sy;
            this.x = this.y;
        else
            xnodes = 1:sx;
            ynodes = longerSideNodes;
            this.y = this.x;
        end
        this.u_t0 = this.u_t0(xnodes, ynodes);
        this.dudt_t0 = this.dudt_t0(xnodes, ynodes);
    end
    
    function [ this, tt, max_v, t, EN, II, v1u, v2u ]= BESolver( this )
        EN = 0;
        II = 0;
        %function [tt, max_v, t, v1u, v2u]=BE2D_v6(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0)
        
        
        h = this.h;
        tau = this.tau;
        beta1 = this.beta1;
        beta2 = this.beta2;
        al = this.alpha;
        estep = this.estep;
        u_t0 = this.u_t0;
        dudt_t0 = this.dudt_t0;

        sy = length(this.y);
        sx = length(this.x);
        order = 2;
        d2x_s = BEUtilities.GetFinDiffMat(7,2,h);
        domainUtils = BEDomainUtils( 1:10, 1:10, order );
        d2x_ss = d2x_s*d2x_s;
        eye7 = eye(7);
        sqD = sqrt(  ( 1 + beta1/tau^2)^2 - 4*beta2/tau^2 );
        X1 =  (1 + beta1/tau^2 + sqD)/( 2 * h^2 / tau^2);
        X2 =  (1 + beta1/tau^2 - sqD)/( 2 * h^2 / tau^2);
        %(beta2/h^2)*(d2x_s-X1*eye7)*(d2x_s-X2*eye7)
        A1 = (-d2x_s+X1*eye7);
        A2 = (-d2x_s+X2*eye7);
        A = h^2*eye7/tau^2 - d2x_s + beta2*d2x_ss/h^2 - beta1*d2x_s/tau^2;
        sA = A(3,1:5);
        sA1 = A1(2,1:3);
        sA2 = A2(2,1:3);
        coeff = (h^2/beta2);
        a11 = A(1,1);

        %Iy = eye(sy);
        W = eigvect(sx);
        IDH = zeros(1,sx);
        for i = 1:sx
           IDH(i) = W(:,i)' * BEUtilities.BandMatMult([(-1) (h^2/beta1+4)  (-1)], W(:,i)', (h^2/beta1+4))';
        end

        max_v(1,1) = max(max(abs(u_t0)));
        max_v(2,1) = max(max(abs(u_t0)));

        [D2x]=BEUtilities.GetFinDiffMat(sx,2,h);
        [D2y]=BEUtilities.GetFinDiffMat(sy,2,h);

        VV = zeros(sx,sy);
        vdah = zeros(sx,sy);
        %tic
        %D2x*u_t0;
        %u_t0*D2y;       
        %toc
        %tic
        %diff2x(u_t0,vdah);
        %diff2y(u_t0,vdah);
        %toc
        b = VV;      b2 = VV;
        v1z = VV; v1mo = VV; v1u = VV;
        v2z = VV; v2mo = VV; v2mt = VV; v2u = VV;        

        [d2vz, d3vz, d4vz] = this.calc_der2d_v2(u_t0,dudt_t0,W,IDH,3,vdah,h,sx,al,beta1,beta2,4);

        v1z = u_t0 + tau*dudt_t0 + (tau^2/2)*d2vz + (tau^3/6)*d3vz + (tau^4/24)*d4vz;

        v2z = v1z;
        v1mo = u_t0; v2mo = u_t0;
        v2mt =  u_t0 - tau*dudt_t0 + (tau^2/2)*d2vz - (tau^3/6)*d3vz + (tau^4/24)*d4vz;
        v1mt = v2mt;
        %figure(1)
        %mesh(x,y,v1z')
        clear('W'); clear('u_t0'); clear('dudt_t0');

        e=1;
        t(1)=0;tt(1)=tau;t(2)=tau;
        k=2;
        %{
        clear('rt');tav = .08;rt(1)=0;
        op=1;
       for o = 2:201
           op=op+1;
          rt(op) = (op-1)*tav;
       end
       %}
    while(t(k)<this.tEnd-1.0e-08)

        b2 = (2*v1z - v1mo)/tau^2;
        VV = domainUtils.XDerivativeZeroBnd ( al* v1z.^2 );
        d2vz = domainUtils.YDerivativeZeroBnd( al* v2z.^2 );
        
        gm = (h^2/beta2)*( 1 + beta1/tau^2);
        d3vz = domainUtils.XDerivativeZeroBnd ( -(beta2/h^2) * domainUtils.YDerivativeZeroBnd(v2z+v1z) - beta1*b2 );
        d4vz = domainUtils.YDerivativeZeroBnd( (-beta2/h^2) * domainUtils.YDerivativeZeroBnd(v2z) + v2z + (beta1/tau^2)*(v2z - 2*v2mo + v2mt) );
        
        b = b2*h^2 + VV + d2vz + d3vz + d4vz;

        f1=0; f2=0;
        for j=1:(sy+1)/2
            if( norm( b(:,j) - b(:,end - j + 1), 2) < 1.00e-10 )
                %v1u(:,j) = diag3solv(sA1,b(:,j));  v1u(:,j) = coeff*diag3solv(sA2,v1u(:,j));
                v1u(:,j) = BEUtilities.PentSolv(a11, sA, b(:,j));
                v1u(:,end - j + 1) = v1u(:,j);%pentsolv(a11,sA,b(:,end - j + 1));
                f1 = f1+1;
            else
                v1u(:,j) = BEUtilities.PentSolv(a11, sA, b(:,j));
                v1u(:,end - j + 1) = BEUtilities.PentSolv(a11, sA, b(:,end - j + 1));
               
                %v1u(:,j) = diag3solv(sA1,b(:,j));  v1u(:,j) = coeff*diag3solv(sA2,v1u(:,j));
                %v1u(:,end - j + 1) = diag3solv(sA1,b(:,end - j + 1));
                %v1u(:,end - j + 1) = coeff*diag3solv(sA2,v1u(:,end - j + 1));
                f2 = f2+1;
            end
        end

        f1=0; f2=0;
        gm = (h^2/beta2)*( 1 + beta1/tau^2);
        b = (2*v2z - v2mo)/tau^2;
        
        %d2vz = al*diff2y((v2z.^2)  , vdah );
        %VV = diff2x (  al* v1z.^2 , vdah );

        d3vz = domainUtils.XDerivativeZeroBnd( ( -beta2/h^2) * domainUtils.XDerivativeZeroBnd(v1u) + (1 + beta1/tau^2)*v1u -( beta2/h^2) * domainUtils.YDerivativeZeroBnd(v2z+v1u ) - beta1*b2 );
        d4vz = domainUtils.YDerivativeZeroBnd(  -beta1 * b );
        b2 = b*h^2  + VV +  d2vz + d3vz + d4vz;
  
        for i=1:(sx+1)/2
            if( norm(b2(i,:) - b2(end - i + 1,:)  , 2) < 1.00e-10 )
                v2u(i,:) = BEUtilities.PentSolv(a11, sA, b2(i,:));
                %v2u(i,:) = diag3solv(sA1,b2(i,:));    v2u(i,:) = coeff*diag3solv(sA2,v2u(i,:));
                v2u(end - i + 1,:) = v2u(i,:);%pentsolv(a11,sA,b2(end - i + 1,:));
                f1 = f1+1;
            else
                v2u(i,:) = BEUtilities.PentSolv(a11, sA, b2(i,:));
                v2u(end - i + 1,:) = BEUtilities.PentSolv(a11, sA, b2(end - i + 1,:));
                
                %v2u(i,:) = diag3solv(sA1,b2(i,:));    v2u(i,:) = coeff*diag3solv(sA2,v2u(i,:));
                %v2u(end - i + 1,:) = diag3solv(sA1,b2(end - i + 1,:));
                %v2u(end - i + 1,:) = coeff*diag3solv(sA2,v2u(end - i + 1,:));
                f2 = f2+1;
            end
        end
                     
        IS_NaN = max(max(isnan(v1u))) + max(max(isnan(v2u)));
        max_v(1,k) = max(max(abs(v1u)));
        max_v(2,k) = max(max(abs(v2u)));
        if(mod(k,estep)==0)
            tt(e)=k*tau;
            %figure(2)
            %mesh(x,y,(vdah + dhb*vz)');
            %title('DH * v')
            %xlabel('x');            ylabel('y');
            %figure(5)
            %mesh(x,y,v1u');
            %mesh(x,y,((v1u - v1mo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            %figure(6)
            %mesh(x,y,v2u');
            %mesh(x,y,((v2u - v2mo)/(2*tau))');
            %title('v_t')
            %xlabel('x');            ylabel('y');
            %colorbar;
            
            str = ['SOL\MM_' num2str(tt(e)) '.mat'];
            save (str,  'v1u');%(1:xstep:end,1:xstep:end);
            str = ['SOL\MM2_' num2str(tt(e)) '.mat'];
            save (str,  'v2u');%(1:xstep:end,1:xstep:end);
            %II(e)=sum(vz)*h;
            %EL(e)=LE(vmo,vz,vu,sdh,sIdh,sdh11,sIdh11,h,tau,sgm);
            %E(e) = EL(e) + NLE(vz,vu,sx,alpha,beta);
            %if(abs(tt(e)-t_end) < 1.00e-010)
            %       tt(end)
            %       t(end)
            %end
            e=e+1;
            
        end
        if( abs( ( this.tEnd-tau + 1.00e-010) - t(k))  < tau )
             if(abs(tt(end)-this.tEnd) > 1.00e-010)
                warning('t_end missed!'); 
                tt(end)
                t(end)      
             else
                 SOL_SAVED = 1
             end
        end
        if(IS_NaN > 0)
            warning('ERROR; BE2D NaN values; Stopping! ');
            t(end) 
            return;
        end
            k=k+1;  
           
            t(k)=(k-1)*tau;
            if((t(k)-floor(t(k)))==0)
                yoyo=t(k)
            end
            v1mo = v1z;
            v1z  = v1u;

            v2mt = v2mo;
            v1mt = v1mo;
            v2mo = v2z;
            v2z  = v2u;
        end
    
        clear('v1mo');  clear('v2mo'); 
        sol_size = size(v1u)
   
    end
    
    function [ name ] = GetName( this )
        name = 'AlternatingMethod';
    end
    
    function [d2vz, d3vz, d4vz, d5vz] = calc_der2d_v2(this, vz,dvz,W,IDH,s_isdh,vdah,h,sx,al,beta1,beta2,order)

        domainUtils = BEDomainUtils( 1:10, 1:10, 2 );
        %==================
        VV = vdah;

        b=(al*vz.*vz + (beta2/h^2)*(-this.dh_gmI(vz,h^2/beta2,domainUtils)));

        deltab = W'* this.dh_gmI(b,0,domainUtils)/beta1; 
        if(s_isdh == 3)  
            for j=1:sx
                VV(j,:) = BEUtilities.TridiagSolv([(-1) IDH(j) (-1)],deltab(j,:));
            end
        else 
        end
        d2vz = W*VV ;
        if(order == 2) d3vz = vdah; d4vz = vdah; d5vz = vdah; return; end

        %==================    

        b=( 2*al*dvz.*vz + (beta2/h^2)*(-this.dh_gmI(dvz,h^2/beta2,domainUtils)));
        deltab = W'* this.dh_gmI(b,0,domainUtils)/beta1;    

        if(s_isdh == 3)  
            for j=1:sx
                VV(j,:) = BEUtilities.TridiagSolv([(-1) IDH(j) (-1)],deltab(j,:));
            end
        else 
        end
        d3vz = W*VV ;
        if(order == 3) d4vz = vdah; d5vz = vdah; return; end

        %==================    

        b=(2*al*(dvz.*dvz + vz.*d2vz) + (beta2/h^2)*(-this.dh_gmI(d2vz,h^2/beta2,domainUtils)));
        deltab = W'* this.dh_gmI(b,0,domainUtils)/beta1;    

        if(s_isdh == 3)  
            for j=1:sx
                VV(j,:) = BEUtilities.TridiagSolv([(-1) IDH(j) (-1)],deltab(j,:));
            end
        else  
        end
        d4vz = W*VV ;
        if(order == 4) d5vz = vdah; return; end

        %==================  

        b=(2*al*(3*dvz.*d2vz + vz.*d3vz) + (beta2/h^2)*(-this.dh_gmI(d3vz, h^2/beta2 ,domainUtils)));
        deltab = W'* this.dh_gmI(b,0,domainUtils)/beta1;    

        if(s_isdh == 3)  
            for j=1:sx
                VV(j,:) = BEUtilities.TridiagSolv([(-1) IDH(j) (-1)],deltab(j,:));
            end
        else  
        end
        d5vz = W*VV ;
        %==================  
    end

    function vdah=dh_gmI(this, M, gm, domainUtils)   %Delta_h - gm*I 
        vdah = (domainUtils.DeltaHZeroBnd(M) - gm*M);
    end
    
    function PlotResults(this, plotBnd, vl1, vl2, t, max_vv)

        if(nargin == 1)
            plotBnd = 0;
        end
         
        viewTypeX = 0;
        viewTypeY = 90;

        Q = 41;
        i = 50;
        if (plotBnd)
            figure(i+1)
            mesh(this.x, this.y(1:Q), vl2(:,1:Q)');
            view( viewTypeX, viewTypeY );
            title('Bottom domain boundary near y=-40');
            xlabel('x');            ylabel('y');
            figure(i+2)
            mesh(this.x(1:Q), this.y, vl2(1:Q,:)');
            view( viewTypeX, viewTypeY );    
            title('Left domain boundary near x=-40');
            xlabel('x');            ylabel('y');

            figure(i+3)
            mesh(this.x, this.y(end-Q:end), vl2(:,end-Q:end)');
            view( viewTypeX, viewTypeY );    
            title('Top domain boundary near y=40');
            xlabel('x');            ylabel('y');    
            figure(i+4)
            mesh(this.x(end-Q:end), this.y, vl2(end-Q:end,:)');
            view( viewTypeX, viewTypeY );
            title('Right domain boundary near x=40');
            xlabel('x');            ylabel('y');
        end

        figure(i+5)
        mesh(this.x,this.y,vl1')
        view( viewTypeX, viewTypeY );
        title('solution 1');
        xlabel('x');            ylabel('y');
        colorbar;
        figure(i+15)
        mesh(this.x,this.y,vl2')
        view( viewTypeX, viewTypeY );
        title('solution 2');
        xlabel('x');            ylabel('y');
        colorbar;

        figure(i+6)
        %hold on;
        plot(t(1:end-1),max_vv(1,:),'k', t(1:end-1),max_vv(2,:),'g', t(1), max_vv(1,1)+0.005, 'k', t(1), max_vv(1,1)-0.005, 'k' );
        %hold off;
        title('Evolution of the maximum');
        xlabel('time "t"');  ylabel('max(u_h)');

    end
    
  end

end

      

        