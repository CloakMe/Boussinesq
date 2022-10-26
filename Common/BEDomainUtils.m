classdef (ConstructOnLoad) BEDomainUtils
    
    properties (SetAccess = private, GetAccess = public)
        X
        Y      
        hx
        hy

        order
        x
        y
        mFiniteDiff
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
      
    function this = BEDomainUtils( x, y, order, der ) 
        
        if(nargin == 3)
            der = 2;
        end
        
        if( mod( order, 2 ) ~= 0 || order<1 )
            error( 'order must be positive even number! ' );
        end
        this.order = order;

        this.hx = x( 2 ) - x( 1 );
        this.hy = y( 2 ) - y( 1 );

        numberOfExtPoints = order/2;
        leftJointX = (-numberOfExtPoints:1:-1)*this.hx + x( 1 );
        rightJointX = (1:numberOfExtPoints)*this.hx + x( end );
        extendedX = [leftJointX x rightJointX];
        
        leftJointY = (-numberOfExtPoints:1:-1)*this.hy + y( 1 );
        rightJointY = (1:numberOfExtPoints)*this.hy + y( end );
        extendedY = [leftJointY y rightJointY];
        [ this.X_xAugDomain, this.Y_xAugDomain ] = this.GetNet(extendedX,y);
        [ this.X_yAugDomain, this.Y_yAugDomain ] = this.GetNet(x,extendedY);
        this.x = x;
        this.y = y;
        [ this.X, this.Y ] = this.GetNet(x,y);
        
        if(der == 2)
            this.mFiniteDiff = BEUtilities.GetFinDiffMatZeroBnd(order+1,order);
        else
            this.mFiniteDiff = BEUtilities.GetFinDiffMat1DerZeroBnd(order+1,order);
        end

%         n = (order/2)-1;
%         for i=1:(order/2)-1
%             pivot = this.mFiniteDiff(this.order/2-i+1, this.order-i+1);
%             for j=1:n
%                 this.mFiniteDiff(j, :) = ...
%                     ( pivot / this.mFiniteDiff(j, this.order-i+1) ) * this.mFiniteDiff(j,:)  -...
%                     this.mFiniteDiff(this.order/2 - i + 1, :);
%             end
%             n = n-1;
%         end
%         this.mFiniteDiff( (order/2)+2:end, :) = flipud( fliplr( this.mFiniteDiff(1:(order/2),:) ) );
        
%       symetric = BEUtilities.GetFinDiffMat(order+1,order);
        boolbreak = 1;
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
    
    function [ zeroMatrix ] = DeltaHZeroBnd( this, M)
        zeroMatrix = this.YDerivativeZeroBnd( M ) + this.YDerivativeZeroBnd(M')';
    end
                                         
    function [ zeroMatrix ] = YDerivativeZeroBnd( this, M ) 

        zeroMatrix = zeros( size( M ) );
        %left boundary elements
        mid = (this.order+2)/2;
        sizeFD = size( this.mFiniteDiff, 2 ) - 1;
        for j = 1:mid-1
            zeroMatrix(:,j) = M(:,1:sizeFD)*this.mFiniteDiff(j,1:end-1)';
        end
        %internal and right boundary elements
        [ zeroMatrix ] = YDerivativeCommonPart( this, M, zeroMatrix, mid, sizeFD );

    end
    
    function [ zeroMatrix ] = XDerivativeZeroBnd( this, M ) 

        zeroMatrix = this.YDerivativeZeroBnd(M')';

    end
    
    function [ zeroMatrix ] = DeltaHEvenFunZeroBnd( this, M)
        zeroMatrix = this.YDerivativeEvenFunZeroBnd( M ) + this.YDerivativeEvenFunZeroBnd(M')';
    end
    
    function [ zeroMatrix ] = YDerivativeEvenFunZeroBnd( this, M ) 

        zeroMatrix = zeros( size( M ) );
        %left boundary elements
        mid = (this.order+2)/2;
        sizeFD = size( this.mFiniteDiff, 2 ) - 1;
        for j = 1:mid-1
            zeroMatrix(:,j) =  M(:,1:mid-1+j)*this.mFiniteDiff(mid,mid+1-j:end)'  + M(:,2:mid+1-j)*this.mFiniteDiff(mid,mid+j:end)';
        end
        
        %internal and right boundary elements
        [ zeroMatrix ] = YDerivativeCommonPart( this, M, zeroMatrix, mid, sizeFD );

    end
    
    function [ zeroMatrix ] = XDerivativeEvenFunZeroBnd( this, M ) 

        zeroMatrix = this.YDerivativeEvenFunZeroBnd(M')';

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
    %internal and right boundary elements
    function [ zeroMatrix ] = YDerivativeCommonPart( this, M, zeroMatrix, mid, sizeFD ) 
       
        for j=mid:size(M,2)-mid+1
            zeroMatrix(:,j) = M(:,j-mid+1:j+mid-1)*this.mFiniteDiff(mid,:)';
        end
        jfd = mid+1;
        for j=size(M,2)-mid+2:size(M,2)
            zeroMatrix(:,j) = M(:,end-sizeFD:end)*this.mFiniteDiff(jfd,1:end)';
            jfd = jfd+1;
        end

    end

  end 

end
