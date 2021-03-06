classdef (ConstructOnLoad) BEDomainUtils
    
    properties (SetAccess = private, GetAccess = protected)
        X
        Y      
        hx
        hy

        order
        x
        y
    end
    
    properties (SetAccess = private, GetAccess = public)
        % this X and Y are net over
        % x augmented domain ( e.g. x(-2), x(-1), x(0), x(1), ... x(end-1), x(end), x(end+1), x(end+2), x(end+3)
        % and normal y, i.e. y(1), y(2), ..., y(end)
        X_xAugDomain 
        Y_xAugDomain
        % this X and Y are net over
        % y augmented domain ( e.g. y(-2), y(-1), y(0), y(1), ... y(end-1), y(end), y(end+1), y(end+2), y(end+3)
        % and normal x, i.e. x(1), x(2), ..., x(end)
        X_yAugDomain
        Y_yAugDomain
    end 
    
  methods( Access = public )
      
    function this = BEDomainUtils( x, y, order ) 
        if( mod( order, 2 ) ~= 0 || order<1 )
            error( 'order must be positive even number! ' );
        end
        this.order = order;
        numberOfExtPoints = order/2;
        this.hx = x( 2 ) - x( 1 );
        this.hy = y( 2 ) - y( 1 );
        
        leftJointX = (-numberOfExtPoints:1:-1)*this.hx + x( 1 );
        rightJointX = (1:numberOfExtPoints)*this.hx + x( end );
        extendedX = [leftJointX x rightJointX];
        
        leftJointY = (-numberOfExtPoints:1:-1)*this.hy + y( 1 );
        rightJointY = (1:numberOfExtPoints)*this.hy + y( end );
        extendedY = [leftJointY y rightJointY];
        
        this.x = x;
        this.y = y;
        [ this.X, this.Y ] = this.GetNet(x,y);
        [ this.X_xAugDomain, this.Y_xAugDomain ] = this.GetNet(extendedX,y);
        [ this.X_yAugDomain, this.Y_yAugDomain ] = this.GetNet(x,extendedY); 
        %this = SetBoundingBox( this, x, y, numberOfExtPoints );

    end
    
    function [ zeroMatrix ] = DeltaH( this,...
                                         M,...
                                         finiteDiff,...
                                         augPntsLeft,...
                                         augPntsRight,...
                                         augPntsTop,...
                                         augPntsBtm ) 

        if( nargin == 3 )
            bndPntsCount = this.order/2;
            augZeroPointsY = zeros( size(M,1), bndPntsCount );
            augZeroPointsX = zeros( bndPntsCount, size(M,2) );
            zeroMatrix = this.YDerivative( M, augZeroPointsY, augZeroPointsY, finiteDiff ) +...
                         this.XDerivative( M, augZeroPointsX, augZeroPointsX, finiteDiff );
        else
            zeroMatrix = this.YDerivative( M, augPntsLeft, augPntsRight, finiteDiff ) +...
                         this.XDerivative( M, augPntsTop, augPntsBtm, finiteDiff );
        end
    end
    
    function [ zeroMatrix ] = XDer(      this,...
                                         M,...
                                         finiteDiff ) 

        augZeroPointsX = zeros( size( this.topX ) );
        zeroMatrix = this.XDerivative( M, augZeroPointsX, augZeroPointsX, finiteDiff );
    end
    
    function [ zeroMatrix ] = YDer(     this,...
                                         M,...
                                         finiteDiff ) 

        augZeroPointsY = zeros( size( this.leftY ) );
        zeroMatrix = this.YDerivative( M, augZeroPointsY, augZeroPointsY, finiteDiff );
    end
    
    function [ zeroMatrix ] = YDerivative( this, M, augLeftPoints, augRightPoints, finiteDiff ) 
        if( size( finiteDiff, 1 ) > 1 )
            error( 'Finite Difference expected as a row vector but received column vector!' );
        end
        [ sxM, syM ] = size( M );
        sizeFD = length( finiteDiff );
        [ sxAugLeftPoints, syAugLeftPoints ] = size( augLeftPoints );
        [ sxAugRightPoints, syAugRightPoints ] = size( augRightPoints );
        
        if( sxM ~= sxAugLeftPoints || sxM ~= sxAugRightPoints )
           error( 'First dimension of augmented matrices is wrong!' ); 
        end        
        if( (sizeFD-1)/2 ~= syAugLeftPoints || (sizeFD-1)/2 ~= syAugRightPoints )
           error( 'Second dimension of augmented matrices is wrong!' ); 
        end

        zeroMatrix = zeros( size( M ) );
        newM = [ augLeftPoints M augRightPoints ];

        for j=1:size(M,2)
            zeroMatrix(:,j) = newM(:,j:j+sizeFD-1)*finiteDiff';
        end

    end
    
    function [ zeroMatrix ] = XDerivative( this, M, augTopPoints, augBtmPoints, finiteDiff ) 
        if( size( finiteDiff, 1 ) > 1 )
            error( 'Finite Difference expected as a row vector but received column vector!' );
        end
        zeroMatrix = this.YDerivative(M',augTopPoints',augBtmPoints',finiteDiff)';

    end
    
    function [X,Y]=GetNet(this,x,y)
        
        sx = length(x);
        sy = length(y);
        ox1 = ones(1,sy);
        X = x'*ox1;
        if(sx~=sy )
            oy1 = ones(1,sx);
            Y = (y'*oy1)';
        else
            Y = (y'*ox1)';
        end
    end
    
  end
  
  methods( Access = protected )
      
    function [ this ] = SetBoundingBox( this, x, y, numberOfExtPoints )
        %================= left ======================
        startAugleftY = y( 1 ) - numberOfExtPoints*this.hy;
        endAugleftY = y( 1 ) - this.hy;
        %x(1)
        %x(end)
        
        yAug = startAugleftY : this.hy : endAugleftY;
        [ this.leftX, this.leftY ] = this.GetNet(x,yAug);
        %=======================================
        %================= right ======================
        startAugrightY = y( end ) + this.hy;
        endAugrightY = y( end ) + numberOfExtPoints*this.hy;
        %x(1)
        %x(end)
        
        yAug = startAugrightY : this.hy : endAugrightY;
        [ this.rightX, this.rightY ] = this.GetNet(x,yAug);
        %=======================================        
        %================== top =====================
        startAugtopX = x( 1 ) - numberOfExtPoints*this.hx;
        endAugtopX = x( 1 ) - this.hx;
        %y(1)
        %y(end)
        
        xAug = startAugtopX : this.hx : endAugtopX;
        [ this.topX, this.topY ] = this.GetNet(xAug,y);
        %=======================================        
        %================== btm =====================
        startAugbtmX = x( end ) + this.hx;
        endAugbtmX = x( end ) + numberOfExtPoints*this.hx;
        %y(1)
        %y(end)
        
        xAug = startAugbtmX : this.hx : endAugbtmX;
        [ this.btmX, this.btmY ] = this.GetNet(xAug,y);
    end
            
    function [X,Y]=GetBndNet(this,x,y, numberOfExtPoints )
        
        sx = length(x);
        sy = length(y);
        
        ox1 = [ ones(1,numberOfExtPoints); zeros(1,sy); ones(1,numberOfExtPoints) ];
        X = x'*ox1;
        if(sx~=sy )
            oy1 = [ ones(1,numberOfExtPoints); zeros(1,sx); ones(1,numberOfExtPoints) ];
            Y = (y'*oy1)';
        else
            Y = (y'*ox1)';
        end
    end
    
  end 

end
