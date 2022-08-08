classdef (ConstructOnLoad) BEEngineTaylorZeroBnd < BEEngineTaylor
   % Class help goes here
  properties ( SetAccess = protected, GetAccess = public)
    domainUtils
  end 
  
  methods
    
    function this = BEEngineTaylorZeroBnd( dscrtParams, eqParams, ic )
        this = this@BEEngineTaylor( dscrtParams, eqParams, ic, 1 );
        Delta_hx = BEUtilities.GetFinDiffMatZeroBnd( this.sx, this.order, this.h );
        Delta_hy = BEUtilities.GetFinDiffMatZeroBnd( this.sy, this.order, this.h );
        
%         n = (this.order/2)-1;
%         for i=1:(this.order/2)-1
%             pivot = Delta_hx(this.order/2-i+1, this.order-i+1);
%             for j=1:n
%                 Delta_hx(j, :) = ...
%                     ( pivot / Delta_hx(j, this.order-i+1) ) * Delta_hx(j,:)  -...
%                     Delta_hx(this.order/2 - i + 1, :);
%             end
%             n = n-1;
%         end
%         Delta_hx( end-(this.order/2)+1:end, :) = flipud( fliplr( Delta_hx(1:(this.order/2),:) ) );
% 
%         n = (this.order/2)-1;
%         for i=1:(this.order/2)-1
%             pivot = Delta_hy(this.order/2-i+1, this.order-i+1);
%             for j=1:n
%                 Delta_hy(j, :) = ...
%                     ( pivot / Delta_hy(j, this.order-i+1) ) * Delta_hy(j,:)  -...
%                     Delta_hy(this.order/2 - i + 1, :);
%             end
%             n = n-1;
%         end        
%         Delta_hy( end-(this.order/2)+1:end, :) = flipud( fliplr( Delta_hy(1:(this.order/2),:) ) );
                
        [ this.eigenFinDiffMat_x, this.D_x ] = eig( Delta_hx' );
        [ this.eigenFinDiffMat_y, this.D_y ] = eig( Delta_hy' );
        this.invEigenFinDiffMat_xT = inv(this.eigenFinDiffMat_x');
        this.invEigenFinDiffMat_y = inv(this.eigenFinDiffMat_y);
        [this.Lambda_X, this.Lambda_Y] = this.SetLambdaNet();
        this.domainUtils = BEDomainUtils( this.x, this.y, this.order );
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
        if(this.beta > 1 || this.beta < 1)
            deltab = this.eigenFinDiffMat_x' * this.domainUtils.DeltaHZeroBnd( nonlinTerm ) * this.eigenFinDiffMat_y;
            X = -deltab ./ ( this.Lambda_X + this.Lambda_Y - this.h^2);
            dnvz = ( this.invEigenFinDiffMat_xT * X * this.invEigenFinDiffMat_y + this.domainUtils.DeltaHZeroBnd( dnU_dtn ) / this.h^2 )/this.beta;
        else
            deltab = this.h^2 * ( this.eigenFinDiffMat_x' * nonlinTerm * this.eigenFinDiffMat_y );
            X = -deltab ./ ( this.Lambda_X + this.Lambda_Y - this.h^2);
            dnvz = ( - nonlinTerm + this.invEigenFinDiffMat_xT * X * this.invEigenFinDiffMat_y + this.domainUtils.DeltaHZeroBnd( dnU_dtn ) / this.h^2 )/this.beta;
        end
        
       
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
        X = -wvt ./ (this.Lambda_X + this.Lambda_Y);
        vec1 = this.h^2 * this.invEigenFinDiffMat_xT * X * this.invEigenFinDiffMat_y;
        
        midv = (vz  + vpo)/2;
        energyTerm = this.beta * ( vec1 + vt ) .* vt  + ...
            (this.beta * midv - this.domainUtils.DeltaHZeroBnd( midv )/this.h^2 ) .* midv + ...
            (this.alpha*this.beta/3 ) * ( vz.^3 + vpo.^3 ); %  
        e = this.GetIntegralOf(energyTerm);
    end
    
    function [lambda_X, lambda_Y] = SetLambdaNet(this)
        
        sx = size(this.D_x,1);
        sy = size(this.D_y,1);
        
        ox1 = ones(1,sy);
        lambda_X = diag(this.D_x)*ox1;
        
        if(sx~=sy )
            oy1 = ones(1,sx);
            lambda_Y = (diag(this.D_y)*oy1)';
        else
            lambda_Y = (diag(this.D_y)*ox1)';
        end
    end
      
  end
  
  properties ( SetAccess = protected, GetAccess = public)
    %Discretization parameters
    eigenFinDiffMat_x
    invEigenFinDiffMat_xT
    D_x
    eigenFinDiffMat_y
    invEigenFinDiffMat_y
    D_y
    Lambda_X
    Lambda_Y
  end 

end