classdef WaveFactory
  methods
    function this = WaveFactory( workspaceName, bndCutSize, nat_42 )
        if( nargin == 2 )
            this = this.WaveFactory1( workspaceName, bndCutSize );
        end
        if( nargin == 1 )
            this = this.WaveFactory2( workspaceName );
        end
        if( nargin == 3 )
            this = this.WaveFactory3( workspaceName );
        end
        
    end
    
    function this = WaveFactory1( this, workspaceName, bndCutSize )
        
        mydir  = pwd;
        idcs   = strfind(mydir,'/');
        if(size(idcs,1) == 0 && size(idcs,2) == 0)
            idcs   = strfind(mydir,'\');
        end
        workspacePath = mydir(1:idcs(end-2));
        workspacePath = strcat(workspacePath,'BEEliptic\Boussinesq2D\SavedWorkspaces\');
        file = fullfile(workspacePath, strcat(workspaceName,'.mat'));
        load( file, '-mat');
        this.x = x;
        this.y = y;
        this.compBox = compBox;
        this.zeroY = zeroY;
        if(length(tauVector)<iterMax && UseExtendedDomain == 1 && size(bigUTimeDerivative,1)~=1)
            this.y_st = y_st2;
            this.y_end = y_end2;
            this.x_end = x_end2;
            this.x_st = x_st2;
        else
            this.y_st = y_st;
            this.y_end = y_end;
            this.x_end = x_end;
            this.x_st = x_st;
        end
        this.h = h;
        this.u_t0 = bigU;
        this.theta = theta(end);
        this.mu = ( c1 + c2 ) / 2.0;
        this.c = c;
        this.beta = bt1/bt2;
        this.beta1 = bt1;
        this.beta2 = bt2;
        this.alpha = al;
        this.order = length( derivative.second ) - 1;
        this.dudt_t0 = this.GetTDer();
        %this.mu = 0;
        bndPtsRem = bndCutSize/this.h;
        this.x = this.x( bndPtsRem+1:end-bndPtsRem );
        this.y = this.y( bndPtsRem+1:end-bndPtsRem );
        this.dudt_t0 = this.dudt_t0( bndPtsRem+1:end-bndPtsRem, bndPtsRem+1:end-bndPtsRem );
        this.u_t0 = this.u_t0( bndPtsRem+1:end-bndPtsRem, bndPtsRem+1:end-bndPtsRem );
    end
    
    function this = WaveFactory2( this, bestFitIc )
        this.mu = 0;
        this.theta = 1;
        run( bestFitIc );
        this.compBox = compBox;
        this.x = x;
        this.y = y;
        this.y_st = y_st;
        this.y_end = y_end;
        this.x_end = x_end;
        this.x_st = x_st;
        this.h = h;
        this.c = c;
        this.beta = beta1/beta2;
        this.beta1 = beta1;
        this.beta2 = beta2;
        this.alpha = al;
        this.order = order;
        if(vc==1)
            this.u_t0 = u_ex2d_mat_vc(X,Y,c,this.beta1);
            if(c==0)
                this.dudt_t0 = 0*this.u_t0;
            else
                tic
                this.dudt_t0 =dudt2d_mat3_vc_v2(X,Y,c,this.beta1);
                toc
            end
        else
            if(vc==0)
            this.u_t0 = u_ex2d_mat_v2(X,Y,this.c,this.beta1);
                if(c==0)
                    this.dudt_t0 = 0*u_t0;
                else
                    tic
                    this.dudt_t0 =dudt2d_mat3_v2(X,Y,c,this.beta1);
                    toc
                end
            end
        end
        if(mod(x(end),h) == 0 && mod(y(end),h) == 0)
            [yo,indx] =min(abs(x));
            [yo,indy] =min(abs(y));
            this.dudt_t0(indx,indy)
            this.dudt_t0(indx,indy) = 0.000;
        end
        IS_dudt_NaN = max(max(isnan(this.dudt_t0)))
        IS_u_NaN = max(max(isnan(this.u_t0)))
        for k = 1:sx
            if( -10^(-11)<x(k) )
                %this.ZeroX = k;
                break;
            end
        end
        for k = 1:sy
            if( -10^(-11)<y(k))
                this.zeroY = k;
                break;
            end
        end
    end
    
    function this = WaveFactory3( this, nat42Ic )
        this.mu = 0;
        this.theta = 1;
        run( nat42Ic );
        this.compBox = compBox;
        this.x = x;
        this.y = y;
        this.y_st = y_st;
        this.y_end = y_end;
        this.x_end = x_end;
        this.x_st = x_st;
        this.h = h;
        this.c = c;
        this.beta = beta1/beta2;
        this.beta1 = beta1;
        this.beta2 = beta2;
        this.alpha = al;
        this.order = order;
        
        this.u_t0 = pola2cart_v5( this.x, this.y, gama1, gama2 );
        if(c==0)
            this.dudt_t0 = 0*this.u_t0;
        else
            tic
            domainUtils = BEDomainUtils( this.x, this.y, this.order );
            this.dudt_t0 = sqrt( -domainUtils.DeltaH( this.u_t0 ) );
            toc
        end
            
        if(mod(x(end),h) == 0 && mod(y(end),h) == 0)
            [yo,indx] =min(abs(x));
            [yo,indy] =min(abs(y));
            this.dudt_t0(indx,indy)
            this.dudt_t0(indx,indy) = 0.000;
        end
        IS_dudt_NaN = max(max(isnan(this.dudt_t0)))
        IS_u_NaN = max(max(isnan(this.u_t0)))
        for k = 1:sx
            if( -10^(-11)<x(k) )
                %this.ZeroX = k;
                break;
            end
        end
        for k = 1:sy
            if( -10^(-11)<y(k))
                this.zeroY = k;
                break;
            end
        end
    end
    
  end
    
    properties ( SetAccess = private, GetAccess = public )
        x
        y
        zeroY
        y_st
        y_end
        x_end
        x_st
        h        
        u_t0
        dudt_t0
        mu
        theta
        alpha
        beta
        beta1
        beta2
        c
        order

        %helpers
        newY
        compBox
    end
    
    
    methods
        function [ clonedWaves, clonedWavesTimeDer ] =...
                CloneWave( this, waveYShift )
            
            colisn = this.zeroY + waveYShift/this.h;

            clonedWaves = [this.u_t0(:,1:colisn) fliplr(this.u_t0(:,1:colisn-1))];
            clonedWavesTimeDer = [this.dudt_t0(:,1:colisn)...
                                  fliplr(this.dudt_t0(:,1:colisn-1))];
        end
        
        function PlotWaves( this, newY, clonedWaves )
            figure(2)
            mesh( this.x, newY, clonedWaves' );
            xlabel('x');    ylabel('y');
            title('two solitons')
            %axis([x_st2 x_end2 y_st2 y_end2 -0.00001 .00001]);
            axis([this.x_st this.x_end newY(1) newY(end) -2 2]);
            colorbar;
            caxis([-0.0001 .0001]);
            view(0,90);
        end
        
        function PlotSingleWave( this, turnOnCaxis )
            figure(2)
            mesh( this.x, this.y, this.u_t0' );
            xlabel('x');    ylabel('y');
            title('soliton')
            %axis([x_st2 x_end2 y_st2 y_end2 -0.00001 .00001]);
            axis([this.x_st this.x_end this.y_st this.y_end -2 2]);
            if( turnOnCaxis ~= 0 )
                colorbar;
                caxis([-0.0001 .0001]);
                view(0,90);
            end
        end
        
        function [ newY ] = GetNewY( this, waveYShift )
            colisn = this.zeroY + waveYShift/this.h;
            newY = -(this.y(colisn)-this.y_st):this.h:(this.y(colisn)-this.y_st);
        end
        
    end
    
    methods( Access = private )
        function dudt_t0 = GetTDer( this )
            
            fdStartPos = this.order/2;
            finiteDiff = BEUtilities.GetFinDiffCoeff( -fdStartPos:fdStartPos, 1)';
            
            domUtils = BEBoundaryUtils( this.x,...
                                        this.y,...
                                        this.order,...
                                        this.beta,...
                                        this.c,...
                                        this.mu,...
                                        this.theta, ...
                                        this.h);
                                    
            %dersBndLeft = domUtils.GetAsymptoticFunctionLeft( 0 );
            %dersBndRight = domUtils.GetDersBndRight( 0 );
            bndPntsCount = this.order/2;
            xNet = domUtils.X_yAugDomain(:,1:bndPntsCount);
            yNet = domUtils.Y_yAugDomain(:,1:bndPntsCount);
            asymptoticFunctionYStart = domUtils.AsymptoticFunction( xNet, yNet, 0, 0 );      
            xNet = domUtils.X_yAugDomain(:,end - bndPntsCount + 1:end);
            yNet = domUtils.Y_yAugDomain(:,end - bndPntsCount + 1:end);
            asymptoticFunctionYEnd = domUtils.AsymptoticFunction( xNet, yNet, 0, 0 );

            dudt_t0 =  -this.c * domUtils.YDerivative( this.u_t0, asymptoticFunctionYStart, asymptoticFunctionYEnd, finiteDiff )/this.h;
            %{
            dudt_t02 =  -this.c * domUtils.YDerivative( this.u_t0, 0*asymptoticFunctionYStart, 0*asymptoticFunctionYEnd, finiteDiff )/this.h;
            figure(1)
            mesh(this.x, this.y(1:7), dudt_t0(:,1:7)');
            xlabel('x');            ylabel('y');
            figure(2)
            mesh(this.x, this.y(1:7), dudt_t02(:,1:7)');
            xlabel('x');            ylabel('y');
            
            fg = 5;
            %}
        end
    end
end
          