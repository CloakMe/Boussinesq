classdef (ConstructOnLoad) BEEngineEnergySaveSoftBnd < BEEngineShrink
   % Class help goes here
  
  methods
    
    function this = BEEngineEnergySaveSoftBnd( dscrtParams, eqParams, ic )
        this = this@BEEngineShrink( dscrtParams, eqParams, ic );
    end
    
    function [this, tt, max_v, t, EN, II, vu, dtv]=BESolver( this )
        %sx2 = sx*sy;
        max_v(1) = max(max(abs(this.u_t0)));
        vu = this.vdah;
        dtv = this.vdah;
        domUtilsEdges = BEDomainUtilsP2Edges( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
        [d2vz, d3vz, d4vz, d5vz, d6vz] = this.Calc_der2d( this.u_t0, this.dudt_t0, 0 );
        
        vu = this.u_t0 + this.tau*this.dudt_t0 + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
            (this.tau^4/24)*d4vz + (this.tau^5/120)*d5vz + (this.tau^6/720)*d6vz; %#ok<*CPROP>
        
        vu = domUtilsEdges.GetExtrapolationOfVzEdges( vu, this.tau );
        %dtv = this.dudt_t0 + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz +...
        %    (this.tau^4/24)*d5vz + (this.tau^5/120)*d6vz;

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
            
            [vu, numOfIter] = this.GetVuPicardi( vz, vmo, t(k) );
            if( numOfIter > 25 )
                return;
            end
                
            %vu = domUtilsEdges.GetExtrapolationOfVzEdges( vu, t(k) );
            IS_dudt_NaN = max(max(isnan(vu)));
            max_v(k) = max(max(abs(vu)));
            if(mod(k,this.estep)==0)
                tt(e)=k*this.tau;
                %{
                %==========================================
                %}
                EN(e) = this.GetEnergy( vz, vu, t( k ) );

                II(e)=sum(sum(vz))*this.h^2;
                
                this.SaveSolutionOnIterStep( tt(e), vu, vz );
                
                e=e+1;
            end
            if(IS_dudt_NaN == 1)
               warning('ERROR; BE2D NaN values; Stopping! ');
               return;
            end
            numAuxPnts = this.order/2;
            [flag, this.maximumOnBnd] = this.CheckForDivergenceOnBnd( vz, numAuxPnts );
            if( flag == 1 )
                warning( 'Algorithm diverges on boundary, shrinking the domain ... ' );
                this = this.UpdateMatricesSizes( this.sx - 2*numAuxPnts, this.sy - 2*numAuxPnts );
                [ vmo, vz, vu ] = this.UpdateMatrices( 1, vmo, vz, vu );
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
        d2vz = this.vdah; d3vz = this.vdah; d4vz = this.vdah; d5vz = this.vdah; d6vz = this.vdah;

        domainUtilsP2edges = BEDomainUtilsP2Edges( this.x,...
                                                 this.y,...
                                                 this.order,...
                                                 this.beta,...
                                                 this.c,...
                                                 this.mu,...
                                                 this.theta );
        domainUtilsP2edges = domainUtilsP2edges.SetDeltaTimeDerivatives( vz, t );
        
        %================== 2nd
        nonlinTerm = this.alpha*this.beta*vz.*vz;
        d2vz = this.GetCurrentDer( vz, nonlinTerm, t, 0, domainUtilsP2edges );
        if( this.order == 2) return; end
 
        %==================  3rd
        nonlinTerm = 2 * this.alpha*this.beta*dvz.*vz;
        d3vz = this.GetCurrentDer( dvz, nonlinTerm, t, 1, domainUtilsP2edges );
        if( this.order == 3) return; end

        %==================  4th
        nonlinTerm = 2*this.alpha*this.beta*(dvz.*dvz + vz.*d2vz);
        d4vz = this.GetCurrentDer( d2vz, nonlinTerm, t, 2, domainUtilsP2edges );
        if( this.order == 4) return; end

        %==================  5th
        nonlinTerm = 2*this.alpha*this.beta*(3*dvz.*d2vz + vz.*d3vz);
        d5vz = this.GetCurrentDer( d3vz, nonlinTerm, t, 3, domainUtilsP2edges );
        if( this.order == 5) d6vz = vdah; return; end
        
        %==================  6th
        nonlinTerm = 2*this.alpha*this.beta*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz);
        d6vz = this.GetCurrentDer( d4vz, nonlinTerm, t, 4, domainUtilsP2edges );
    
    end
    
    function [dnvz] = GetCurrentDer( this,...
                                     timeDerVz,...
                                     nonlinTerm,...
                                     t,...
                                     timeDerOrd,...
                                     domainUtilsP2spec)
        if( this.order > 4 )
            error( 'Not yet implemented for order = 6 (or 8) !' );
        end      
        VV = this.vdah;
        fd2ndDer = this.GetFd2ndDer();
        mid = ( this.order/2 + 1 );
        deltab = this.eigenFinDiffMat'*(...
            domainUtilsP2spec.DeltaH( nonlinTerm, fd2ndDer ) + ...
            (this.beta-1) * domainUtilsP2spec.DeltaTimeDerevative( timeDerVz, fd2ndDer, timeDerOrd ) );    
        deltab = ( deltab )/this.beta;
        for j=1:this.sx
            diag = [ -fd2ndDer(1:mid-1) this.IminusDHdiag(j) -fd2ndDer(mid+1:end) ];
            %diag = [(1/12) (-16/12) this.IminusDHdiag(j) (-16/12) 1/12];
            if( this.order == 2 )      
                VV(j,:) = BEUtilities.TridiagSolv( diag, deltab(j,:) );
            end
            if( this.order == 4 )   
                VV(j,:) = BEUtilities.PentSolv( this.IminusDHdiag(j), diag, deltab(j,:));
            end
        end
        
        deltav = ( domainUtilsP2spec.DeltaTimeDerevative( timeDerVz, fd2ndDer, timeDerOrd ) )/this.beta;
        dnvz = ( this.eigenFinDiffMat*VV + deltav/this.h^2 );
    end
       
    function [ vu, numOfIter ] = GetVuPicardi( this, vz, vmo, t )
        eps = 10^(-13);
        
        %spparms('bandden',0.5);
    %===============================================
        domainUtilsP2edges = ...
            BEDomainUtilsP2Edges( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
        domainUtilsP2edges = domainUtilsP2edges.SetDeltaTimeDerivatives( vz, t );
            
        niv0 = vz;
        %===============================================
        nonlinTerm = this.alpha*this.beta*( vz.^2 + vz.*vmo + vmo.^2)/3;
        d2vz =  this.GetCurrentDer( vz,  nonlinTerm, t, 2, domainUtilsP2edges );
        vu = 2*vz - vmo + this.tau^2 * d2vz;
        vu = domainUtilsP2edges.GetExtrapolationOfVzEdges( vu, this.tau );
        %===============================================
        mmax=max( max( abs( niv0 ) ) );
        numOfIter=0;
        
        while( max( max( abs( niv0 - vu ) ) ) > eps * mmax )
            niv0=vu;
            mmax = max( abs( vu ) );
            if( numOfIter > 25 )
                warning('ERROR; too much iterations in SIT! ');
                return;
            end
            numOfIter=numOfIter+1; %5
        %===============================================
            nonlinTerm =  this.alpha*this.beta*( niv0.^2 + niv0.*vmo + vmo.^2 )/3;
            d2vz =  this.GetCurrentDer( vz,  nonlinTerm, t, 2, domainUtilsP2edges );
            vu = 2*vz - vmo + this.tau^2 * d2vz; 
            vu = domainUtilsP2edges.GetExtrapolationOfVzEdges( vu, this.tau );
        %===============================================
            %max(abs(niv0-niv));
         %   niv(1:20)
        end
    end
  end

end