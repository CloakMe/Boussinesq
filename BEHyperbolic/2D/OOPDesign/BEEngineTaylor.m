classdef (ConstructOnLoad) BEEngineTaylor < BEEngine
   % Class help goes here
  
  methods
    
    function this = BEEngineTaylor( dscrtParams, eqParams, ic )
        this = this@BEEngine( dscrtParams, eqParams, ic );
    end
    
    function [ this, tt, max_v, t, EN, II, vu, dtv ]= BESolver( this )
        %sx2 = sx*sy;
        max_v(1) = max(max(abs(this.u_t0)));
        s_idh = 5;        
        vu = this.vdah;

        %[d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(this.u_t0,this.dudt_t0,eigenFinDiffMat,IminusDHdiag,s_idh,vdah, this.h, this.sx,this.alpha,this.beta,this.order);
        [d2vz, d3vz, d4vz, d5vz, d6vz] = this.Calc_der2d( this.u_t0, this.dudt_t0, 0 );

        vu = this.u_t0 + this.tau*this.dudt_t0 + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
            (this.tau^4/24)*d4vz + (this.tau^5/120)*d5vz + (this.tau^6/720)*d6vz; %#ok<*CPROP>
        dtv = this.dudt_t0 + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz +...
            (this.tau^4/24)*d5vz + (this.tau^5/120)*d6vz;

        %figure(2)
        %mesh(x,y,(this.tau^2/2)*d2vz')
        %title('(this.tau^2/2)*d^2v/dt^2 , t=this.tau')
        %xlabel('x');            ylabel('y');
        clear('this.dudt_t0'); clear('d2vz');

        e=1; tt(e) = 0;
        t(1)=0;t(2)=this.tau;
        k=2;
        vz = vu; vmo = this.u_t0;
        clear('this.u_t0'); clear('v2');
        EN(1)=0;II(1)=0;
        while( t(k) < this.tEnd )
            %[d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vu,dtv,eigenFinDiffMat,IminusDHdiag,s_idh,vdah, this.h, this.sx,this.alpha,this.beta,this.order);
            [d2vz, d3vz, d4vz, d5vz, d6vz] = this.Calc_der2d( vu, dtv, t(k) );
            vu = vu + this.tau*dtv + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
                (this.tau^4/24)*d4vz+ (this.tau^5/120)*d5vz + (this.tau^6/720)*d6vz;
            dtv = dtv + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz +...
                (this.tau^4/24)*d5vz + (this.tau^5/120)*d6vz;    

            IS_dudt_NaN = max(max(isnan(vu)));
            max_v(k) = max(max(abs(vu)));
            if(mod(k,this.estep)==0)
                tt(e)=k*this.tau;
                
                EN(e) = this.GetEnergy( vz, vu, t( k ) );
                %EN(e) = Eg2_vc_2d(vmo, vz, vu, this.minusDHdiag, this.eigenFinDiffMat,this.h,this.tau,...
                %        this.alpha,this.beta, this.sx, this.sy,s_idh, this.vdah,this.vdah,this.vdah);
                this.SaveSolutionOnIterStep( tt(e), vu, vz );
                II(e)=sum(sum(vz))*this.h^2;
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
                 else
                     SOL_SAVED = 1
                 end
            end
                k=k+1;  

                t(k)=(k-1)*this.tau;
                if((t(k)-floor(t(k)))==0)
                    yoyo=t(k)
                end
            vmo = vz; vz = vu; 
        end
        
        clear('W');    clear('vmo'); 
        sol_size = size(vu) 
    end

    function [d2vz, d3vz, d4vz, d5vz, d6vz] = Calc_der2d( this, vz,dvz, t )
        if( this.order > 6 )
            error( 'Max derivative order is 6th!' ); 
        end
        %[ augDhDomainOnly ] = this.GetDhAugDomainOnly( t );
        %[ augDtDomainOnly ] = this.GetDtAugDomainOnly( t );
        %[ augPowDomainOnly ] = this.GetPowAugDomainOnly( t );
        d2vz = this.vdah; d3vz = this.vdah; d4vz = this.vdah; d5vz = this.vdah; d6vz = this.vdah;

        %================== 2nd
        nonlinTerm = (this.alpha*this.beta*vz.*vz +...
                     (this.beta-1)*vz);
        d2vz = this.GetCurrentDer( nonlinTerm, vz, t, 0 );%, augDhDomainOnly(:,:,1), augDtDomainOnly(:,:,1) );
        if( this.order == 2) return; end
 
        %==================  3rd   
        nonlinTerm =(2*this.alpha*this.beta*dvz.*vz +...
                    (this.beta-1)*dvz);
        d3vz = this.GetCurrentDer( nonlinTerm, dvz, t, 1);%, augDhDomainOnly(:,:,2), augDtDomainOnly(:,:,2) );
        if( this.order == 3) return; end

        %==================  4th  
        nonlinTerm =( 2*this.alpha*this.beta*(dvz.*dvz + vz.*d2vz) +...
                    (this.beta-1)*d2vz );
        d4vz = this.GetCurrentDer( nonlinTerm, d2vz, t, 2);%, augDhDomainOnly(:,:,3), augDtDomainOnly(:,:,3) );
        if( this.order == 4) return; end

        %==================  5th
        nonlinTerm =( 2*this.alpha*this.beta*(3*dvz.*d2vz + vz.*d3vz) +...
                    (this.beta-1)*d3vz );
        d5vz = this.GetCurrentDer( nonlinTerm, d3vz, t, 3);%, augDhDomainOnly(:,:,4), augDtDomainOnly(:,:,4) );
        if( this.order == 5)  d6vz = vdah; return; end
        
        %==================  6th
        nonlinTerm =( 2*this.alpha*this.beta*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) +...
                    (this.beta-1)*d4vz );
        d6vz = this.GetCurrentDer( nonlinTerm, d4vz, t, 4);% augDhDomainOnly(:,:,5), augDtDomainOnly(:,:,5) );
    
    end
    
    function [dnvz] = GetCurrentDer( this,...
                                     nonlinTerm,...
                                     highestDer,...
                                     t, ...
                                     order)
        if( this.order > 6 )
            error( 'Not yet implemented for order = 8!' );
        end
        boundaryUtils = BEBoundaryUtils( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
                
        VV = this.vdah;
        VV2 = this.vdah;
        fd2ndDer = this.GetFd2ndDer();
        mid = ( this.order/2 + 1 );
        %{
        bndPntsCount = this.order/2;
        xNet = boundaryUtils.X_yAugDomain(:,1:bndPntsCount);
        yNet = boundaryUtils.Y_yAugDomain(:,1:bndPntsCount);
        asymptoticFunctionYStart = boundaryUtils.AsymptoticFunction( xNet, yNet, t, order );      
        xNet = boundaryUtils.X_yAugDomain(:,end - bndPntsCount + 1:end);
        yNet = boundaryUtils.Y_yAugDomain(:,end - bndPntsCount + 1:end);
        asymptoticFunctionYEnd = boundaryUtils.AsymptoticFunction( xNet, yNet, t, order );

        dhighestDerdt2 =  this.c^2 * boundaryUtils.YDerivative( VV, asymptoticFunctionYStart, asymptoticFunctionYEnd, fd2ndDer );
        %}        
        deltab = this.eigenFinDiffMat'* ( boundaryUtils.DeltaH( nonlinTerm, fd2ndDer )+...
            boundaryUtils.DeltaXH( this.beta-1, fd2ndDer, t, order ) );    

        %deltab2 = this.eigenFinDiffMat'* ( boundaryUtils.DeltaH( nonlinTerm, fd2ndDer ));  
        
        for j=1:this.sx
            diag = [ -fd2ndDer(1:mid-1) this.IminusDHdiag(j) -fd2ndDer(mid+1:end) ];
            %diag = [(1/12) (-16/12) this.IminusDHdiag(j) (-16/12) 1/12];
            if( this.order == 2 )      
                VV(j,:) = BEUtilities.TridiagSolv( diag, deltab(j,:) );
            end
            if( this.order == 4 )   
                VV(j,:) = BEUtilities.PentSolv( this.IminusDHdiag(j), diag, deltab(j,:));
                %VV2(j,:) = BEUtilities.SevenSolv( this.IminusDHdiag(j), diag, deltab2(j,:));
            end
            if( this.order == 6 )
                VV(j,:) = BEUtilities.SevenSolv( this.IminusDHdiag(j), diag, deltab(j,:));
                %VV2(j,:) = BEUtilities.SevenSolv( this.IminusDHdiag(j), diag, deltab2(j,:));
            end
        end
        
        deltav = ( boundaryUtils.DeltaH( highestDer, fd2ndDer ) )/this.h^2;%/this.beta;
        myBoundary = boundaryUtils.DeltaXH( 1, fd2ndDer, t, order )/this.h^2;%
        dnvz = ( this.eigenFinDiffMat*VV + deltav + myBoundary)/this.beta;
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
       
  end

end