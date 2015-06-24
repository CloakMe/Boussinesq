classdef WaveFactory
    methods
        function this = WaveFactory( workspaceName, bndPtsRem )
            workspaceName = strcat( 'D:\workspace\Matlab\Boussinesq\BEEliptic\Boussinesq2D\SavedWorkspaces\', workspaceName, '.mat');
            load( workspaceName, '-mat');
            
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
            this.beta1 = bt1;
            this.beta2 = bt2;
            this.alpha = al;
            this.order = length( derivative.second ) - 1;
            
            this.dudt_t0 = this.GetTDer();
            %this.mu = 0;
            
            this.x = this.x( bndPtsRem+1:end-bndPtsRem );
            this.y = this.y( bndPtsRem+1:end-bndPtsRem );
            this.dudt_t0 = this.dudt_t0( bndPtsRem+1:end-bndPtsRem, bndPtsRem+1:end-bndPtsRem );
            this.u_t0 = this.u_t0( bndPtsRem+1:end-bndPtsRem, bndPtsRem+1:end-bndPtsRem );
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
            fd = BEUtilities.GetFinDiffCoeff( -fdStartPos:fdStartPos, 1);
            
            domUtils = BEDomainUtilsP2( this.x,...
                                        this.y,...
                                        this.order,...
                                        this.beta1/this.beta2,...
                                        this.c,...
                                        this.mu,...
                                        this.theta );
                                    
             dersBndLeft = domUtils.GetDersBndLeft( 0 );
             dersBndRight = domUtils.GetDersBndRight( 0 );
                 
             dudt_t0 = -this.c * domUtils.YDerivative(...
                                             this.u_t0,...
                                             dersBndLeft(:,:,1),...
                                             dersBndRight(:,:,1),...
                                             fd' )/this.h;
        end
    end
end
          