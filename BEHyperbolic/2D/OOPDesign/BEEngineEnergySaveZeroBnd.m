classdef (ConstructOnLoad) BEEngineEnergySaveZeroBnd < BEEngine
   % Class help goes here
  
  methods
    
    function this = BEEngineEnergySaveZeroBnd( dscrtParams, eqParams, ic )
        dscrtParams.order = 2;
        this = this@BEEngine( dscrtParams, eqParams, ic );
    end
    
    function [this, tt, max_v, t, EN, II, vu, dtv]=BESolver( this )
        %sx2 = sx*sy;
        max_v(1) = max(max(abs(this.u_t0)));
        vu = this.vdah;
        dtv = this.vdah;
        % BEDomainUtilsP2Edges( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
        % domUtilsEdges = BEDomainUtils( this.x, this.y, this.order); 
        this.order = 6;
        [d2vz, d3vz, d4vz, d5vz, d6vz] = this.Calc_der2d( this.u_t0, this.dudt_t0, 0 );
        
        vu = this.u_t0 + this.tau*this.dudt_t0 + (this.tau^2/2)*d2vz + (this.tau^3/6)*d3vz +...
            (this.tau^4/24)*d4vz + (this.tau^5/120)*d5vz + (this.tau^6/720)*d6vz; %#ok<*CPROP>
        this.order = 2;
        %dtv = this.dudt_t0 + this.tau*d2vz + (this.tau^2/2)*d3vz + (this.tau^3/6)*d4vz +...
        %    (this.tau^4/24)*d5vz + (this.tau^5/120)*d6vz;
        
        EN_1 = this.GetEnergy( this.u_t0, vu, 0 );
        II_1 = this.GetIntegralOf( this.u_t0 );
        
        %figure(2)
        %mesh(x,y,(this.tau^2/2)*d2vz')
        %title('(this.tau^2/2)*d^2v/dt^2 , t=this.tau')
        %xlabel('x');            ylabel('y');
        clear('this.dudt_t0'); clear('d2vz');

        e=2;
        t = zeros(1,this.tEnd/this.tau);
        tt = zeros(1,this.tEnd/this.tau);
        t(1)=0;t(2)=this.tau;
        k=2;
        vz = vu; vmo = this.u_t0;
        clear('this.u_t0'); clear('v2');
        EN=zeros(1,this.tEnd/this.tau);
        II=zeros(1,this.tEnd/this.tau);
        EN(1) = EN_1;
        II(1) = II_1;
        while( t(k) < this.tEnd )
            
            [vu, numOfIter] = this.GetVuPicardi( vz, vmo, t(k) );
            if( numOfIter > 25 )
                return;
            end
                
            %vu = domUtilsEdges.GetExtrapolationOfVzEdges( vu, t(k) );
            IS_dudt_NaN = max(max(isnan(vu)));
            max_v(k) = max(max(abs(vu)));
            EN(k) = this.GetEnergy( vz, vu, t( k ) );
            II(k)= this.GetIntegralOf( vz );
            if(mod(k,this.estep)==0)
                tt(e)=k*this.tau;
                %{
                mesh(this.x, this.y, vu')
                title('vu , t=this.tau')
                xlabel('x');            ylabel('y');
                %==========================================
                %}               
                %this.SaveSolutionOnIterStep( tt(e), vu );                
                e=e+1;
            end
            if(IS_dudt_NaN == 1)
               warning('ERROR; BE2D NaN values; Stopping! ');
               return;
            end

            if( abs( ( this.tEnd - this.tau + 1.00e-010 ) - t(k) )  < this.tau )
                 if( abs( tt(end) - this.tEnd ) > 1.00e-010 )
                    warning( 'End time missed!' ); 
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

    function [d2vz, d3vz, d4vz, d5vz, d6vz] = Calc_der2d( this, vz,dvz, t )
        if( this.order > 6 )
            error( 'Max derivative order is 6th!' ); 
        end
        d2vz = this.vdah; d3vz = this.vdah; d4vz = this.vdah; d5vz = this.vdah; d6vz = this.vdah;
        % BEDomainUtilsP2Edges( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
        domainUtils = BEDomainUtils( this.x,...
                                                 this.y,...
                                                 this.order ); %BEDomainUtilsP2Edges
        %domainUtils = domainUtils.SetDeltaTimeDerivatives( vz, t );
        
        %================== 2nd
        nonlinTerm = this.alpha*this.beta*vz.*vz + ( this.beta - 1 )*vz;
        d2vz = this.GetCurrentDer( vz, nonlinTerm, t, 0, domainUtils );
        if( this.order == 2) return; end
 
        %==================  3rd
        nonlinTerm = 2 * this.alpha*this.beta*dvz.*vz + ( this.beta - 1 )*dvz;
        d3vz = this.GetCurrentDer( dvz, nonlinTerm, t, 1, domainUtils );
        if( this.order == 3) return; end

        %==================  4th
        nonlinTerm = 2*this.alpha*this.beta*(dvz.*dvz + vz.*d2vz) + ( this.beta - 1 )*d2vz;
        d4vz = this.GetCurrentDer( d2vz, nonlinTerm, t, 2, domainUtils );
        if( this.order == 4) return; end

        %==================  5th
        nonlinTerm = 2*this.alpha*this.beta*(3*dvz.*d2vz + vz.*d3vz) + ( this.beta - 1 )*d3vz;
        d5vz = this.GetCurrentDer( d3vz, nonlinTerm, t, 3, domainUtils );
        if( this.order == 5) d6vz = vdah; return; end
        
        %==================  6th
        nonlinTerm = 2*this.alpha*this.beta*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) + ( this.beta - 1 )*d4vz;
        d6vz = this.GetCurrentDer( d4vz, nonlinTerm, t, 4, domainUtils );
    
    end
    
    function [dnvz] = GetCurrentDer( this,...
                                     timeDerVz,...
                                     nonlinTerm,...
                                     t,...
                                     timeDerOrd,...
                                     domainUtils)
        if( this.order > 6 )
            error( 'Not yet implemented for order = 8 !' );
        end      
        
        
        VV = this.vdah;
        fd2ndDer = this.GetFd2ndDer();
        mid = ( this.order/2 + 1 );
        deltab = this.eigenFinDiffMat'*(...
            domainUtils.DeltaH( nonlinTerm, fd2ndDer ) ); %DeltaTimeDerevative 

        for j=1:this.sx
            diag = [ -fd2ndDer(1:mid-1) this.IminusDHdiag(j) -fd2ndDer(mid+1:end) ];
            %diag = [(1/12) (-16/12) this.IminusDHdiag(j) (-16/12) 1/12];
            if( this.order == 2 )
                VV(j,:) = BEUtilities.TridiagSolv( diag, deltab(j,:) );
            elseif( this.order == 4 )
                VV(j,:) = BEUtilities.PentSolv( this.IminusDHdiag(j), diag, deltab(j,:));
            elseif( this.order == 6 )
                VV(j,:) = BEUtilities.SevenSolv( this.IminusDHdiag(j), diag, deltab(j,:));
            end
        end
        
        deltav = ( domainUtils.DeltaH( timeDerVz, fd2ndDer ) ); % DeltaTimeDerevative
        dnvz = ( this.eigenFinDiffMat*VV + deltav/this.h^2 )/this.beta;
    end
       
    function [ vu, numOfIter ] = GetVuPicardi( this, vz, vmo, t )
        eps = 10^(-13);
        
        %spparms('bandden',0.5);
    %===============================================
    % BEDomainUtilsP2Edges( this.x, this.y, this.order, this.beta, this.c, this.mu, this.theta );
        domainUtils = ...
            BEDomainUtils( this.x, this.y, this.order ); %BEDomainUtilsP2Edges
        %sdomainUtils = domainUtils.SetDeltaTimeDerivatives( vz, t );
            
        niv0 = vz;
        %===============================================
        nonlinTerm = this.alpha*this.beta*( vz.^2 + vz.*vmo + vmo.^2)/3 +...
            ( this.beta - 1 ) * ( vz + vmo )/2;
        d2vz =  this.GetCurrentDer( vz ,  nonlinTerm, t, 2, domainUtils );
        vu = 2*vz - vmo + this.tau^2 * d2vz;
        
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
            nonlinTerm =  this.alpha*this.beta*( niv0.^2 + niv0.*vmo + vmo.^2 )/3 +...
                ( this.beta - 1 ) * ( niv0 + vmo )/2;
            d2vz =  this.GetCurrentDer( vz,  nonlinTerm, t, 2, domainUtils );
            vu = 2*vz - vmo + this.tau^2 * d2vz; 
        %===============================================
            %max(abs(niv0-niv));
         %   niv(1:20)
        end
    end
  end

end