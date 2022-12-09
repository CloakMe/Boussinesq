classdef (ConstructOnLoad) BEEngineTaylorLyche < BEEngine
   % Class help goes here
  
  methods
    
    function this = BEEngineTaylorLyche( dscrtParams, eqParams, ic, flag )
        if(nargin == 3)
            flag = 0;
        end
        this = this@BEEngine( dscrtParams, eqParams, ic, flag );
        this.mBoundaryUtils = BEBoundaryUtils( this.x, this.y, this.order, this.alpha, this.beta, this.c, this.mu, this.theta, this.h );
        
        Delta_hx = BEUtilities.GetFinDiffMat( this.sx, this.order, this.h );
        Delta_hy = BEUtilities.GetFinDiffMat( this.sy, this.order, this.h );

        [ this.eigenFinDiffMat_x, this.D_x ] = eig( Delta_hx' );
        [ this.eigenFinDiffMat_y, this.D_y ] = eig( Delta_hy' );
        this.invEigenFinDiffMat_xT = inv(this.eigenFinDiffMat_x');
        this.invEigenFinDiffMat_y = inv(this.eigenFinDiffMat_y);
        [this.Lambda_X, this.Lambda_Y] = this.SetLambdaNet();         
    end
    
    function [ this, tt, max_v, t, EN, II, vu, dtv ]= BESolver( this )
        %sx2 = sx*sy;
        max_v(1) = max(max(abs(this.u_t0)));
        s_idh = 5;        
        vu = this.vdah;

        %[d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(this.u_t0,this.dudt_t0,eigenFinDiffMat,IminusDHdiag,s_idh,vdah, this.h, this.sx,this.alpha,this.beta,this.order);
        [d2vz, d3vz, d4vz, d5vz, d6vz, d7vz] = this.Calc_der2d( this.u_t0, this.dudt_t0, 0 );
        zerofy = [ this.order >= 3, this.order >= 5 ];
        vu = this.u_t0 + this.tau*this.dudt_t0 + (this.tau^2/2)*d2vz + zerofy(1) * (this.tau^3/6)*d3vz +...
            (this.tau^4/24)*d4vz + zerofy(2) * (this.tau^5/120)*d5vz + (this.tau^6/720)*d6vz; %#ok<*CPROP>
        dtv = this.dudt_t0 + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz +...
            (this.tau^4/24)*d5vz + (this.tau^5/120)*d6vz + (this.tau^6/720)*d7vz;
        
        EN_1 = this.GetEnergy( this.u_t0, vu, 0 );
        II_1 = this.GetIntegralOf( this.u_t0 );
        this.SaveSolutionOnIterStep( 0, this.u_t0, this.dudt_t0 ); 
        %figure(2)
        %mesh(x,y,(this.tau^2/2)*d2vz')
        %title('(this.tau^2/2)*d^2v/dt^2 , t=this.tau')
        %xlabel('x');            ylabel('y');
        clear('this.dudt_t0'); clear('d2vz');        
        e=2;
        Tsize = int8(this.tEnd/this.tau);
        t = zeros(1,Tsize);
        tt = zeros(1,Tsize);
        t(1)=0;t(2)=this.tau;
        k=2;
        vz = vu; vmo = this.u_t0;
        clear('this.u_t0'); clear('v2');
        EN = zeros(1,Tsize);
        II=zeros(1,Tsize);        
        EN(1) = EN_1;
        II(1) = II_1;
        while( t(k) < this.tEnd )
            %[d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vu,dtv,eigenFinDiffMat,IminusDHdiag,s_idh,vdah, this.h, this.sx,this.alpha,this.beta,this.order);
            [d2vz, d3vz, d4vz, d5vz, d6vz, d7vz] = this.Calc_der2d( vu, dtv, t(k) );
            vu = vu + this.tau*dtv + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
                (this.tau^4/24)*d4vz + (this.tau^5/120)*d5vz + (this.tau^6/720)*d6vz;
            dtv = dtv + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz +...
                (this.tau^4/24)*d5vz + (this.tau^5/120)*d6vz + (this.tau^6/720)*d7vz;    

            IS_dudt_NaN = max(max(isnan(vu)));
            max_v(k) = max(max(abs(vu)));
            EN(k) = this.GetEnergy( vz, vu, t( k ) );
            II(k)= this.GetIntegralOf( vz );
            if(mod(k,this.estep)==0)
                tt(e)=k*this.tau;
                if( mod(k*this.tau, 1) == 0 )
                    this.SaveSolutionOnIterStep( tt(e), vu, dtv ); 
                end 
                e=e+1;
            end
            if(IS_dudt_NaN == 1)
               warning('ERROR; BE2D NaN values; Stopping! ');
               return;
            end
            if( abs( ( this.tEnd - this.tau + 1.00e-010 ) - t(k) )  < this.tau )
                 if( abs( tt(end) - this.tEnd ) > 1.00e-010 )
                    warning( 'tEnd missed!' ); 
                    tt(end)
                    t(end)
                 end
            end
            k=k+1;

            t(k)=(k-1)*this.tau;
            if((t(k)-floor(t(k)))==0)
                fprintf('iteration time: %.4f\n', t(k));
            end
            vmo = vz; vz = vu; 
        end
        II(k)= this.GetIntegralOf( vu );
        clear('W');    clear('vmo'); 
        sol_size = size(vu) 
    end
    
    function [ name ] = GetName( this )
        name = 'TaylorLyche';
    end
    
    function [d2vz, d3vz, d4vz, d5vz, d6vz, d7vz] = Calc_der2d( this, vz,dvz, t )
        if( this.order > 6 )
            error( 'Max derivative order is 6th!' ); 
        end
        %[ augDhDomainOnly ] = this.GetDhAugDomainOnly( t );
        %[ augDtDomainOnly ] = this.GetDtAugDomainOnly( t );
        %[ augPowDomainOnly ] = this.GetPowAugDomainOnly( t );
        d2vz = this.vdah; d3vz = this.vdah; d4vz = this.vdah; d5vz = this.vdah; d6vz = this.vdah; d7vz = this.vdah;
        
        %================== 2nd
        nonlinTerm = (this.alpha*this.beta*vz.*vz +...
                     (this.beta-1)*vz);
        d2vz = this.GetCurrentDer( nonlinTerm, vz, t, 0 );%, augDhDomainOnly(:,:,1), augDtDomainOnly(:,:,1) );
        if( this.order + 1 == 2) return; end
        %{
        fd2ndDer = this.GetFd2ndDer();
        boundaryUtils = BEBoundaryUtils( this.x, this.y, this.order, this.alpha, this.beta, this.c, this.mu, this.theta, this.h );
         Q=4;
        f1 =  ( this.beta*d2vz - boundaryUtils.DeltaH( vz, fd2ndDer )/this.h^2);
        figure(1)
        mesh(this.x(1+Q:end-Q), this.y(1+Q:1+Q+10), f1(1+Q:end-Q,1+Q:1+Q+10)');
        xlabel('x');            ylabel('y');
        
        f2 =  ( d2vz - boundaryUtils.DeltaH( vz, fd2ndDer )/this.h^2);
        figure(2)
        mesh(this.x(1+Q:end-Q), this.y(1+Q:1+Q+10), f2(1+Q:end-Q,1+Q:1+Q+10)');
        xlabel('x');            ylabel('y');
        
        figure(3)
        mesh(this.x(1+Q:end-Q), this.y(1+Q:1+Q+10), vz(1+Q:end-Q,1+Q:1+Q+10)');
        xlabel('x');            ylabel('y');
        
        fg = 5;
        %}
        %==================  3rd   
        nonlinTerm =(2*this.alpha*this.beta*dvz.*vz +...
                    (this.beta-1)*dvz);
        d3vz = this.GetCurrentDer( nonlinTerm, dvz, t, 1);%, augDhDomainOnly(:,:,2), augDtDomainOnly(:,:,2) );
        if( this.order + 1 == 3) return; end

        %==================  4th  
        nonlinTerm =( 2*this.alpha*this.beta*(dvz.*dvz + vz.*d2vz) +...
                    (this.beta-1)*d2vz );
        d4vz = this.GetCurrentDer( nonlinTerm, d2vz, t, 2);%, augDhDomainOnly(:,:,3), augDtDomainOnly(:,:,3) );
        if( this.order + 1 == 4) return; end

        %==================  5th
        nonlinTerm =( 2*this.alpha*this.beta*(3*dvz.*d2vz + vz.*d3vz) +...
                    (this.beta-1)*d3vz );
        d5vz = this.GetCurrentDer( nonlinTerm, d3vz, t, 3);%, augDhDomainOnly(:,:,4), augDtDomainOnly(:,:,4) );
        if( this.order + 1 == 5) return; end
        
        %==================  6th
        nonlinTerm =( 2*this.alpha*this.beta*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) +...
                    (this.beta-1)*d4vz );
        d6vz = this.GetCurrentDer( nonlinTerm, d4vz, t, 4);% augDhDomainOnly(:,:,5), augDtDomainOnly(:,:,5) );
        if( this.order + 1 == 6) return; end
        
         %==================  7th
        nonlinTerm =( 2*this.alpha*this.beta*(10*d2vz.*d3vz + 5*dvz.*d4vz + vz.*d5vz) +...
                    (this.beta-1)*d5vz );
        d7vz = this.GetCurrentDer( nonlinTerm, d5vz, t, 5);% augDhDomainOnly(:,:,5), augDtDomainOnly(:,:,5) );
    end
    
    function [dnvz] = GetCurrentDer( this,...
                                     nonlinTerm,...
                                     dnU_dtn,...
                                     t, ...
                                     order)
        if( this.order > 6 )
            error( 'Not yet implemented for order = 8!' );
        end

        fd2ndDer = this.GetFd2ndDer();
        
        b1 = this.mBoundaryUtils.AssymptoticFuncSquareOutside( fd2ndDer, t, order ) + ...
        	this.mBoundaryUtils.AssymptoticFuncOutside( this.beta-1, fd2ndDer, t, order );
        
        deltab = this.eigenFinDiffMat_x' * ( this.mBoundaryUtils.DeltaH( nonlinTerm, fd2ndDer ) + b1 ) * this.eigenFinDiffMat_y;
        X = -deltab ./ ( this.Lambda_X + this.Lambda_Y - this.h^2);
        
        myBoundary = this.mBoundaryUtils.AssymptoticFuncOutside( 1, fd2ndDer, t, order )/this.h^2;%+ myBoundary
        
        dnvz =  ( this.invEigenFinDiffMat_xT * X * this.invEigenFinDiffMat_y + this.mBoundaryUtils.DeltaH( dnU_dtn, fd2ndDer )/this.h^2 + myBoundary)/ this.beta;
        
        %max1 = max( max( b1 ) )
        %max2 = max( max( b2 ) )
       
%         for j=1:this.sx
%             diag = [ -fd2ndDer(1:mid-1) this.IminusDHdiag(j) -fd2ndDer(mid+1:end) ];
%             %diag = [(1/12) (-16/12) this.IminusDHdiag(j) (-16/12) 1/12];
%             if( this.order == 2 )      
%                 VV(j,:) = BEUtilities.TridiagSolv( diag, deltab(j,:) );
%             elseif( this.order == 4 )   
%                 VV(j,:) = BEUtilities.PentSolv( this.IminusDHdiag(j), diag, deltab(j,:));
%             elseif( this.order == 6 )
%                 VV(j,:) = BEUtilities.SevenSolv( this.IminusDHdiag(j), diag, deltab(j,:));
%             end
%         end
%         
%         deltav = ( this.mBoundaryUtils.DeltaH( dnU_dtn, fd2ndDer ) )/this.h^2;%/this.beta;
% 
%         myBoundary = this.mBoundaryUtils.AssymptoticFuncOutside( 1, fd2ndDer, t, order )/this.h^2;%+ myBoundary
%         dnvz = ( this.eigenFinDiffMat*VV + deltav + myBoundary)/this.beta;
        %{
        Q = 21;
        dnvz2 =  ( this.eigenFinDiffMat*VV2 + deltav )/this.beta;
        figure(1)
        mesh(this.x, this.y(1:Q), dnvz(:,1:Q)');
        xlabel('x');            ylabel('y');
        figure(2)
        mesh(this.x, this.y(1:Q), dnvz2(:,1:Q)');
        xlabel('x');            ylabel('y');
        
        fg = 5;
        %}
    end
    
    function [e] = GetEnergy( this, vz, vpo, t )
        vt = (vpo - vz)/this.tau;
        
        wvt = this.eigenFinDiffMat_x' * vt * this.eigenFinDiffMat_y;        
        X = -wvt ./ (this.Lambda_X + this.Lambda_Y);
        vec1 = this.h^2 * this.invEigenFinDiffMat_xT * X * this.invEigenFinDiffMat_y;
        
        fd2ndDer = this.GetFd2ndDer();
        energyTerm = this.beta * ( vec1 + vt ) .* vt  + ...
            (this.beta * vz - this.mBoundaryUtils.DeltaH( vz, fd2ndDer)/this.h^2 ) .* vz + ...
            (2*this.alpha*this.beta/3 ) * ( vz.^3 ); %  
        e = this.GetIntegralOf(energyTerm);
    end
      
  end
  
  properties ( SetAccess = protected, GetAccess = public)
    mBoundaryUtils
  end 

end