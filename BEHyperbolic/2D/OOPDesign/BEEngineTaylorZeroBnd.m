classdef (ConstructOnLoad) BEEngineTaylorZeroBnd < BEEngineTaylor
   % Class help goes here
  
  
  methods
    
    function this = BEEngineTaylorZeroBnd( dscrtParams, eqParams, ic )
        this = this@BEEngineTaylor( dscrtParams, eqParams, ic, 1 );
        this.Delta_hx = BEUtilities.GetFinDiffMatZeroBnd( this.sx, this.order, this.h );
        this.Delta_hy = BEUtilities.GetFinDiffMatZeroBnd( this.sy, this.order, this.h );

        [ this.eigenFinDiffMat_x, this.w_x ] = eig( -this.Delta_hx );
        [ this.eigenFinDiffMat_y, this.w_y ] = eig( -this.Delta_hy );
        [this.Lambda_X, this.Lambda_Y] = this.SetLambdaNet();
    end
    
    function [ this, tt, max_v, t, EN, II, vu, dtv ]= BESolver( this )
        [ this, tt, max_v, t, EN, II, vu, dtv ] = BESolver@BEEngineTaylor(this);
    end
    
    function [ name ] = GetName( this )
        name = 'TaylorZeroBnd';
    end
        
    function [dnvz] = GetCurrentDer( this,...
                                     nonlinTerm,...
                                     dnU_dtn,...
                                     t, ...
                                     order)
        if( this.order > 6 )
            error( 'Not yet implemented for order = 8!' );
        end
      
        deltab = this.eigenFinDiffMat_x' * ( this.Delta_hx * nonlinTerm + nonlinTerm * this.Delta_hy ) * this.eigenFinDiffMat_y;    
        X = deltab ./ (this.Lambda_X + this.Lambda_Y + this.h^2);
        dnvz = ( this.eigenFinDiffMat_x * X * this.eigenFinDiffMat_y' + this.Delta_hx / this.h^2 * dnU_dtn + dnU_dtn * this.Delta_hy / this.h^2 )/this.beta;
        %{
        Q = 21;
        dnvz2 =  ( deltab )/this.beta;
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

        idhv = this.beta * (vz+vpo) - ( this.Delta_hx * ( vz+vpo ) + ( vz+vpo ) * this.Delta_hy )/this.h^2;
        X = wvt ./ (this.Lambda_X + this.Lambda_Y);
        vec1 = this.eigenFinDiffMat_x * X * this.eigenFinDiffMat_y';  %/h^2
        sigma = 0;  
        IDhvt = this.beta*vt - ( this.Delta_hx * vt + vt * this.Delta_hy )/this.h^2;
        Le= this.beta * ( vec1 + vt ) .* vt  +...
            this.tau^2 * ( sigma - 1/4 ) *  IDhvt.*vt + idhv .* ( vz+vpo ) / 4 ;

        NLe = ( this.alpha*this.beta/3 ) * ( vz.^3  + vpo.^3 );

        energyTerm = Le + NLe;
        e = this.GetIntegralOf(energyTerm);
    end
    
    function [lambda_X, lambda_Y] = SetLambdaNet(this)
        
        sx = size(this.Delta_hx,1);
        sy = size(this.Delta_hy,1);
        
        ox1 = ones(1,sy);
        lambda_X = diag(this.w_x)*ox1;
        
        if(sx~=sy )
            oy1 = ones(1,sx);
            lambda_Y = (diag(this.w_y)*oy1)';
        else
            lambda_Y = (diag(this.w_y)*ox1)';
        end
    end
      
  end
  
  properties ( SetAccess = protected, GetAccess = public)
    %Discretization parameters
    Delta_hx
    Delta_hy    
    eigenFinDiffMat_x
    w_x
    eigenFinDiffMat_y
    w_y
    Lambda_X
    Lambda_Y
  end 

end