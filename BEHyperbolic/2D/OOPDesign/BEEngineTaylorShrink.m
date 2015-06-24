classdef (ConstructOnLoad) BEEngineTaylorShrink < BEEngineShrink
   % Class help goes here
  
  methods
    
    function this = BEEngineTaylorShrink( dscrtParams, eqParams, ic )
        this = this@BEEngineShrink( dscrtParams, eqParams, ic );
    end
    
    function [this, tt, max_v, t, EN, II, vu, dtv]=BESolver( this )
        %sx2 = sx*sy;
        max_v(1) = max(max(abs(this.u_t0)));
        s_idh = 5;        
        vu = this.vdah;

        %[d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(this.u_t0,this.dudt_t0,eigenFinDiffMat,IminusDHdiag,s_idh,vdah, this.h, this.sx,this.alpha,this.beta,this.order);
        [d2vz, d3vz, d4vz, d5vz, d6vz, this] = this.Calc_der2d( this.u_t0, this.dudt_t0, 0 );
        
        [ this.u_t0, this.dudt_t0 ] = this.UpdateMatrices( this.order/2, this.u_t0, this.dudt_t0 );

        if( this.order == 2 )
            vu = this.u_t0 + this.tau*this.dudt_t0 + (this.tau^2/2)*d2vz; %#ok<*CPROP>
            dtv = this.dudt_t0 + this.tau*d2vz;
        end
        
        if( this.order == 4 )
            %[ d2vz, d3vz ] = this.UpdateMatrices( this.order/2 - 1, d2vz, d3vz );
            vu = this.u_t0 + this.tau*this.dudt_t0 + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
                (this.tau^4/24)*d4vz; %#ok<*CPROP>
            dtv = this.dudt_t0 + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz;
        end
               
        clear('this.dudt_t0'); clear('d2vz');

        e=1; tt(e) = 0;
        t(1)=0;t(2)=this.tau;
        k=2;
        vz = vu; vmo = this.u_t0;
        clear('this.u_t0'); clear('v2');
        EN(1)=0;II(1)=0;
        while(t(k)<this.tEnd)
            %[d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vu,dtv,eigenFinDiffMat,IminusDHdiag,s_idh,vdah, this.h, this.sx,this.alpha,this.beta,this.order);
            [d2vz, d3vz, d4vz, d5vz, d6vz, this] = this.Calc_der2d( vu, dtv, t(k) );
            
            [ vu, dtv, vmo, vz ] = this.UpdateMatrices( this.order/2, vu, dtv, vmo, vz );
            [ this.u_t0, this.dudt_t0 ] = this.UpdateMatrices( this.order/2, this.u_t0, this.dudt_t0 );

            if( this.order == 2 )
                vu = vu + this.tau*dtv + (this.tau^2/2)*d2vz; %#ok<*CPROP>
                dtv = this.dudt_t0 + this.tau*d2vz;
            end

            if( this.order == 4 )
                %[ d2vz, d3vz ] = this.UpdateMatrices( this.order/2 - 1, d2vz, d3vz );
                vu = vu + this.tau*dtv + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
                    (this.tau^4/24)*d4vz; %#ok<*CPROP>
                dtv = this.dudt_t0 + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz;
            end
            
            figure(4)
            mesh( this.x, this.y, (  vu  )' );
            title('vu')
            xlabel('x');            ylabel('y');
            IS_dudt_NaN = max(max(isnan(vu)));
            max_v(k) = max(max(abs(vu)));
            if(mod(k,this.estep)==0)
                tt(e)=k*this.tau;                
                
                %figure(2)
                %mesh(this.x,this.y,d4vz'); %(this.tau^2/2)*
                %title('(this.tau^2/2)*d^2v/dt^2 , t=this.tau')
                %xlabel('x');            ylabel('y');
                EN(e) = this.GetEnergy( vz, vu, t( k ) );
                
                II(e)=sum(sum(vz))*this.h^2;
                
                this.SaveSolutionOnIterStep( tt(e), vu, vz );
                
                e=e+1;
            end
            if(IS_dudt_NaN == 1)
               warning('ERROR; BE2D NaN values; Stopping! ');
               return;
            end
            if( abs( ( this.tEnd-this.tau + 1.00e-010) - t(k))  < this.tau )
                 if(abs(tt(end)-this.tEnd) > 1.00e-010)
                    warning('tEnd missed!'); 
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
    
    function [d2vz, d3vz, d4vz, d5vz, d6vz, this] = Calc_der2d( this, vz,dvz, t )
        if( this.order > 6 )
            error( 'Max derivative order is 6th!' ); 
        end
        
        d2vz = this.vdah; d3vz = this.vdah; d4vz = this.vdah; d5vz = this.vdah; d6vz = this.vdah;

        %================== 2nd
        nonlinTerm = (this.alpha*this.beta*vz.*vz + (this.beta-1)*vz);
        [ d2vz, ty ] = this.GetCurrentDer( nonlinTerm, vz, t );
        if( this.order == 2) 
            this = ty;
            return; 
        end
 
        %==================  3rd   
        nonlinTerm = (2*this.alpha*this.beta*dvz.*vz + (this.beta-1)*dvz);
        [ d3vz, this ] = this.GetCurrentDer( nonlinTerm, dvz, t );
        if( this.order == 3) return; end

        %==================  4th  
        numAuxPnts = this.order / 2;
        vz = vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        dvz = dvz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        nonlinTerm = ( 2*this.alpha*this.beta*(dvz.*dvz + vz.*d2vz) + (this.beta-1)*d2vz );
        
        [ d4vz, gg ] = this.GetCurrentDer( nonlinTerm, d2vz, t );
        if( this.order == 4) 
            this = gg; 
            d2vz = d2vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
            d3vz = d3vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
            return; 
        end

        %==================  5th
        nonlinTerm =( 2*this.alpha*this.beta*(3*dvz.*d2vz + vz.*d3vz) + (this.beta-1)*d3vz );
        [ d5vz, this ] = this.GetCurrentDer( nonlinTerm, d3vz, t );
        if( this.order == 5) return; end
        
        %==================  6th
        vz = vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        dvz = dvz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        d2vz = d2vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        d3vz = d3vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        nonlinTerm =( 2*this.alpha*this.beta*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) +...
                    (this.beta-1)*d4vz );
        d6vz = this.GetCurrentDer( nonlinTerm, d4vz, t );
        
        d2vz = d2vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        d3vz = d3vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        d4vz = d4vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        d5vz = d5vz( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
    end
    
    function [dnvz, this] = GetCurrentDer( this, nonlinTerm, highestDer, t )
        if( this.order > 4 )  
            error( 'Not yet implemented for order = 6 (or 8) !' );
        end
        domainUtilsP2 = BEDomainUtilsP2( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
        
        numAuxPnts = this.order / 2;
        this = this.UpdateMatricesSizes( this.sx - 2*numAuxPnts, this.sy - 2*numAuxPnts );
                
        fd2ndDer = this.GetFd2ndDer();
        mid = ( this.order/2 + 1 );
        deltaHalfB = domainUtilsP2.DeltaH( nonlinTerm, fd2ndDer );
        
        deltab = deltaHalfB( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        deltab = this.eigenFinDiffMat'* deltab;
        %deltab = ( deltab )/this.beta;
        
        VV = this.vdah;
        for j=1:this.sx
            diag = [ -fd2ndDer(1:mid-1) this.IminusDHdiag(j) -fd2ndDer(mid+1:end) ];
            %diag = [(1/12) (-16/12) this.IminusDHdiag(j) (-16/12) 1/12];
            if( this.order == 2 )      
                VV(j,:) = BEUtilities.TridiagSolv( diag, deltab(j,:) );
            end
            if( this.order == 4 )
                VV(j,:) = BEUtilities.PentSolv( this.IminusDHdiag(j), diag, deltab(j,:) );
            end
        end
        deltaHalfv = domainUtilsP2.DeltaH( highestDer, fd2ndDer );
        deltav = deltaHalfv( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        %deltav = ( deltav )/this.beta;
        dnvz = ( this.eigenFinDiffMat*VV + deltav/this.h^2 );
    end
    
  end

end