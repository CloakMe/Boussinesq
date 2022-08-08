classdef (ConstructOnLoad) BEDomainUtilsP2Edges < BEDomainUtilsP2
  properties (SetAccess = private, GetAccess = public)
        % left right top btm are with respect to matrix not domain!
        leftIndecesX
        leftIndecesY
        rightIndecesX
        rightIndecesY
        topIndecesX
        topIndecesY
        btmIndecesX
        btmIndecesY
        %leftDeltaTimeDerivatives
        %rightDeltaTimeDerivatives
        %topDeltaTimeDerivatives
        %btmDeltaTimeDerivatives
        
        topLeftIndecesX
        topLeftIndecesY
        topRightIndecesX
        topRightIndecesY
        btmLeftIndecesX
        btmLeftIndecesY
        btmRightIndecesX
        btmRightIndecesY
        topLeftDeltaTimeDerivatives
        topRightDeltaTimeDerivatives
        btmLeftDeltaTimeDerivatives
        btmRightDeltaTimeDerivatives
        
        topLeftSqIndecesX
        topLeftSqIndecesY
        topRightSqIndecesX
        topRightSqIndecesY
        btmLeftSqIndecesX
        btmLeftSqIndecesY
        btmRightSqIndecesX
        btmRightSqIndecesY
  end
  
  methods
      
    function this = BEDomainUtilsP2Edges( x, y, order, beta, c, mu, theta )
        this = this@BEDomainUtilsP2( x, y, order, beta, c, mu, theta );
        numberOfBndPoints = order/2;
        this = SetIndecesBoundingBox( this, x, y, numberOfBndPoints );
    end
    
    function [ this ] = SetDeltaTimeDerivatives( this, v, t ) %( this, X, Y, v, t, timeDerOrd )
        this.topLeftDeltaTimeDerivatives =...
            this.DeltaTimeDersBnd( this.X( this.topLeftIndecesX, this.topLeftIndecesY ),...
                                   this.Y( this.topLeftIndecesX, this.topLeftIndecesY ),...
                                   v(      this.topLeftIndecesX, this.topLeftIndecesY ),...
                                   t,...
                                   this.order );
        this.topRightDeltaTimeDerivatives =...
            this.DeltaTimeDersBnd( this.X( this.topRightIndecesX, this.topRightIndecesY ),...
                                   this.Y( this.topRightIndecesX, this.topRightIndecesY ),...
                                   v(      this.topRightIndecesX, this.topRightIndecesY ),...
                                   t,...
                                   this.order );
        this.btmLeftDeltaTimeDerivatives =...
            this.DeltaTimeDersBnd( this.X( this.btmLeftIndecesX, this.btmLeftIndecesY ),...
                                   this.Y( this.btmLeftIndecesX, this.btmLeftIndecesY ),...
                                   v(      this.btmLeftIndecesX, this.btmLeftIndecesY ),...
                                   t,...
                                   this.order );
        this.btmRightDeltaTimeDerivatives =...
            this.DeltaTimeDersBnd( this.X( this.btmRightIndecesX, this.btmRightIndecesY ),...
                                   this.Y( this.btmRightIndecesX, this.btmRightIndecesY ),...
                                   v(      this.btmRightIndecesX, this.btmRightIndecesY ),...
                                   t,...
                                   this.order );
    end
        
    function [ zeroMatrix ] = YDer( this, M, finiteDiff )
        if( size( finiteDiff, 1 ) > 1 )
            error( 'Finite Difference expected as a row vector but received column vector!' );
        end
        [ sxM, syM ] = size( M );
        
        zeroMatrix = zeros( size( M ) );
        
        numOfBndPoints = this.order/2;
        for j=numOfBndPoints+1:syM-numOfBndPoints
            zeroMatrix(:,j) = M(:,j-numOfBndPoints:j+numOfBndPoints)*finiteDiff';
        end

    end
    
    function [ zeroMatrix ] = XDer( this, M, finiteDiff )
        [ zeroMatrix ] = YDer( this, M', finiteDiff )';
    end
    
    function [ yder ] = YDeriv( this, M, finiteDiff, timeDerOrd  )
        yder = YDer( this, M, finiteDiff );
        
        yder( this.leftIndecesX, this.leftIndecesY ) =...
            1/(this.c^2-1)*this.XDer( M( this.leftIndecesX, this.leftIndecesY ), finiteDiff );
        yder( this.rightIndecesX, this.rightIndecesY ) =...
            1/(this.c^2-1)*this.XDer( M( this.rightIndecesX, this.rightIndecesY ), finiteDiff );
        
        
        
        numOfBndPoints = this.order/2;
        yder(1:numOfBndPoints,1:numOfBndPoints) = this.topLeftDeltaTimeDerivatives(:,:,timeDerOrd+1);
        yder(end-numOfBndPoints+1:end,1:numOfBndPoints) = this.btmLeftDeltaTimeDerivatives(:,:,timeDerOrd+1);
        yder(1:numOfBndPoints,end-numOfBndPoints+1:end) = this.topRightDeltaTimeDerivatives(:,:,timeDerOrd+1);
        yder(end-numOfBndPoints+1:end,end-numOfBndPoints+1:end) =...
            this.btmRightDeltaTimeDerivatives(:,:,timeDerOrd+1);
        
    end
    
    function [ xder ] = XDeriv( this, M, finiteDiff, timeDerOrd )
        xder = XDer( this, M, finiteDiff );
        
        xder( this.topIndecesX, this.topIndecesY ) =...
            (this.c^2-1)*this.YDer( M(this.topIndecesX, this.topIndecesY), finiteDiff );
        xder( this.btmIndecesX, this.btmIndecesY ) =...
            (this.c^2-1)*this.YDer( M(this.btmIndecesX, this.btmIndecesY), finiteDiff );
        
    end
    
    function [ deltaTimeder ] = DeltaTimeDerevative( this, M, finiteDiff, timeDerOrd  )
        deltaTimeder = YDeriv( this, M, finiteDiff, timeDerOrd ) + ...
                       XDeriv( this, M, finiteDiff, timeDerOrd );        
    end
        
    % === time derivatives
    function deltaTimeDerivatives = DeltaTimeDersBnd( this, X, Y, v, t, timeDerOrd )
        
        if( nargin == 4 )
            timeDerOrd = this.order;
        end
        
        if( max(size(v) ~= size(X)) || max( size(v) ~= size(Y) ) )
            error( 'Function and domain size missmatch!' );
        end
        
        deltaTimeDerivatives = zeros( size( X, 1 ), size( X, 2 ), this.order+1 );
        %{
            figure(9)
          mesh( X(:,1), Y(1,:), v' )
        figure(10)
            mesh( X(:,1), Y(1,:), A' )
        figure(11)
            mesh( X(:,1), Y(1,:), mu' )
        %}
        c = this.c;
        bt = this.beta;
        theta = this.theta;
        omcs = 1 - c^2;
                
        A = omcs*X.^2 - (Y-c*t).^2;
        B = omcs*X.^2 + (Y-c*t).^2;
        bndFun = theta * ( A ./ B.^2 );
        mu = v ./ bndFun;

        %fun = theta * mu .* (  omcs * X.^2 - (Y - c * t).^2 ) ./ ((1-c^2) * X.^2 + (Y - c * t).^2).^2;
        
        deltaTimeDerivatives(:,:,1) = 6 * theta * mu .* (-1 + omcs) .* (c ^ 4 * t ^ 4 - 4 * Y * c ^ 3 * t ^ 3 +...
            6 * Y .^ 2 * c ^ 2 * t ^ 2 - 6 * omcs * X .^ 2 * c ^ 2 * t ^ 2 - 4 * Y .^ 3 * c * t +...
            12 * omcs * X .^ 2 .* Y * c * t + omcs ^ 2 * X .^ 4 - 6 * omcs * X .^ 2 .* Y .^ 2 + Y .^ 4) ./...
            (c ^ 2 * t ^ 2 - 2 * Y * c * t + omcs * X .^ 2 + Y .^ 2).^ 4;
        
        if( timeDerOrd == 0 )
            return;
        end
        
        deltaTimeDerivatives(:,:,2) =  -24 * theta * mu .* c .* (c * t - Y) * (-1 + omcs) .* (5 * omcs ^ 2 * X .^ 4 +...
            20 * omcs * X .^ 2 .* Y * c * t - 10 * omcs * X .^ 2 .* Y .^ 2 - 10 * omcs * X .^ 2 * c ^ 2 * t ^ 2 +...
            c ^ 4 * t ^ 4 - 4 * Y .^ 3 * c * t + 6 * Y .^ 2 * c ^ 2 * t ^ 2 - 4 * Y * c ^ 3 * t ^ 3 + Y .^ 4) ./...
            (c ^ 2 * t ^ 2 - 2 * Y * c * t + omcs * X .^ 2 + Y .^ 2).^ 5;
                
        if( this.order == 1 || timeDerOrd == 1 )
            return;
        end

        deltaTimeDerivatives(:,:,3) = 120 * theta * mu .* c ^ 2 * (-1 + omcs) .* (-omcs * X .^ 2 + Y .^ 2 - 2 * Y * c * t +...
            c ^ 2 * t ^ 2) .* (c ^ 4 * t ^ 4 - 4 * Y * c ^ 3 * t ^ 3 + 6 * Y .^ 2 * c ^ 2 * t ^ 2 -... 
            14 * omcs * X .^ 2 * c ^ 2 * t ^ 2 + 28 * omcs * X .^ 2 .* Y * c * t - 4 * Y .^ 3 * c * t -...
            14 * omcs * X .^ 2 .* Y .^ 2 + Y .^ 4 + omcs ^ 2 * X .^ 4) ./ (c ^ 2 * t ^ 2 - 2 * Y * c * t +...
            omcs * X .^ 2 + Y .^ 2).^ 6;
                
        if( this.order == 2 || timeDerOrd == 2 )
            return;
        end

        deltaTimeDerivatives(:,:,4) = -720 * c ^ 3 * theta * mu .* (c * t - Y) .* (-1 + omcs) .* (-7 * omcs ^ 3 * X .^ 6 -...
            70 * Y .* c * omcs ^ 2 .* X .^ 4 * t + 35 * omcs ^ 2 * X .^ 4 .* Y .^ 2 +...
            35 * c ^ 2 * omcs ^ 2 * X .^ 4 * t ^ 2 + 84 * Y .* c ^ 3 * omcs .* X .^ 2 * t ^ 3 -...
            126 * Y .^ 2 .* c ^ 2 * omcs .* X .^ 2 * t ^ 2 + 84 * Y .^ 3 .* c * omcs .* X .^ 2 * t -...
            21 * Y .^ 4 .* omcs .* X .^ 2 - 21 * c ^ 4 * omcs * X .^ 2 * t ^ 4 + c ^ 6 * t ^ 6 -...
            6 * Y * c ^ 5 * t ^ 5 + Y .^ 6 - 20 * Y .^ 3 * c ^ 3 * t ^ 3 + 15 * Y .^ 2 * c ^ 4 * t ^ 4 +...
            15 * Y .^ 4 * c ^ 2 * t ^ 2 - 6 * Y .^ 5 * c * t) ./ (c ^ 2 * t ^ 2 - 2 * Y * c * t +...
            omcs * X .^ 2 + Y .^ 2).^ 7;
        
        if(this.order == 3 || timeDerOrd == 3 )
            return;
        end

        deltaTimeDerivatives(:,:,5) = 5040 * theta * mu .* c ^ 4 * (-1 + omcs) .* (-8 * Y .^ 7 * c * t -...
            56 * Y .^ 3 * c ^ 5 * t ^ 5 + 70 * Y .^ 4 * c ^ 4 * t ^ 4 +...
            28 * Y .^ 2 * c ^ 6 * t ^ 6 + 28 * Y .^ 6 * c ^ 2 * t ^ 2 - 8 * Y * c ^ 7 * t ^ 7 -...
            56 * Y .^ 5 * c ^ 3 * t ^ 3 - 28 * c ^ 6 * omcs * X .^ 2 * t ^ 6 +...
            70 * c ^ 4 * omcs ^ 2 * X .^ 4 * t ^ 4 - 28 * c ^ 2 * omcs ^ 3 * X .^ 6 * t ^ 2 + Y .^ 8 +...
            c ^ 8 * t ^ 8 + 70 * omcs ^ 2 * X .^ 4 .* Y .^ 4 - 28 * omcs ^ 3 * X .^ 6 .* Y .^ 2 -...
            28 * omcs * X .^ 2 .* Y .^ 6 + 560 * c ^ 3 * Y .^ 3 .* omcs .* X .^ 2 * t ^ 3 -...
            420 * c ^ 4 * Y .^ 2 .* omcs .* X .^ 2 * t ^ 4 + 168 * c ^ 5 * Y .* omcs .* X .^ 2 * t ^ 5 -...
            280 * c ^ 3 * Y .* omcs ^ 2 .* X .^ 4 * t ^ 3 - 420 * c ^ 2 * Y .^ 4 .* omcs .* X .^ 2 * t ^ 2 +...
            420 * c ^ 2 * omcs ^ 2 * X .^ 4 .* Y .^ 2 * t ^ 2 + 168 * Y .^ 5 .* c * omcs .* X .^ 2 * t -...
            280 * Y .^ 3 .* c * omcs ^ 2 .* X .^ 4 * t + 56 * Y .* c * omcs ^ 3 .* X .^ 6 * t + omcs ^ 4 * X .^ 8) ./...
            (c ^ 2 * t ^ 2 - 2 * Y * c * t + omcs * X .^ 2 + Y .^ 2).^ 8;
        if(this.order == 4 || timeDerOrd == 4 )
            return;
        end
        
        if(this.order > 4)
            error( 'Time derivatives for the boundary function are up to 4th order!  ' );
        end
    end
    
    function testFunMid = GetTimeDersBndMid( this, t, derOrd )
        if( nargin == 2 )
            testFunMid = DersTimeBnd( this, this.X, this.Y, t );
        else
            testFunMid = DersTimeBnd( this, this.X, this.Y, t, derOrd );
        end
    end
    
    function testFunLeft = GetTimeDersBndLeft( this, t, derOrd )
        if( nargin == 2 )
            testFunLeft = DersTimeBnd( this, this.leftX, this.leftY, t );
        else
            testFunLeft = DersTimeBnd( this, this.leftX, this.leftY, t, derOrd );
        end
    end
    
    function testFunRight = GetTimeDersBndRight( this, t, derOrd )
        if( nargin == 2 )
            testFunRight = DersTimeBnd( this, this.rightX, this.rightY, t );
        else
            testFunRight = DersTimeBnd( this, this.rightX, this.rightY, t, derOrd );
        end        
    end
    
    function testFunTop = GetTimeDersBndTop( this, t, derOrd )
        if( nargin == 2 )
            testFunTop = DersTimeBnd( this, this.topX, this.topY, t ); 
        else
            testFunTop = DersTimeBnd( this, this.topX, this.topY, t, derOrd );
        end        
    end
    
    function testFunBtm = GetTimeDersBndBtm( this, t, derOrd )
        if( nargin == 2 )
            testFunBtm = DersTimeBnd( this, this.btmX, this.btmY, t );
        else
            testFunBtm = DersTimeBnd( this, this.btmX, this.btmY, t, derOrd );
        end        
    end
    
    function [ vz ] = GetExtrapolationOfVzEdges( this, vz, t )
        
        
        [ topLeftSqX, topLeftSqY ]=this.GetNet( this.x( this.topLeftSqIndecesX ), this.y( this.topLeftSqIndecesY ) );
        [ topRightSqX, topRightSqY ]=this.GetNet( this.x( this.topRightSqIndecesX ), this.y( this.topRightSqIndecesY ) );
        [ btmLeftSqX, btmLeftSqY ]=this.GetNet( this.x( this.btmLeftSqIndecesX ), this.y( this.btmLeftSqIndecesY ) );
        [ btmRightSqX, btmRightSqY ]=this.GetNet( this.x( this.btmRightSqIndecesX ), this.y( this.btmRightSqIndecesY ) );
        
        [ topLeftEdge ] = this.ExtrapolateEdge( vz( this.topLeftSqIndecesX, this.topLeftSqIndecesY ),...
                                                t,...
                                                topLeftSqX,...
                                                topLeftSqY,...
                                                this.topLeftIndecesX,...
                                                this.topLeftIndecesY );
        [ topRightEdge ] = this.ExtrapolateEdge( vz( this.topRightSqIndecesX, this.topRightSqIndecesY ),...
                                                 t,...
                                                 topRightSqX,...
                                                 topRightSqY,...
                                                 this.topRightIndecesX,...
                                                 this.topRightIndecesY);
        [ btmLeftEdge ] = this.ExtrapolateEdge( vz( this.btmLeftSqIndecesX, this.btmLeftSqIndecesY ),...
                                                t,...
                                                btmLeftSqX,...
                                                btmLeftSqY,...
                                                this.btmLeftIndecesX,...
                                                this.btmLeftIndecesY);
        [ btmRightEdge ] = this.ExtrapolateEdge( vz( this.btmRightSqIndecesX, this.btmRightSqIndecesY ),...
                                                 t,...
                                                 btmRightSqX,...
                                                 btmRightSqY,...
                                                 this.btmRightIndecesX,...
                                                 this.btmRightIndecesY );
          
        vz( this.topLeftIndecesX , this.topLeftIndecesY ) = topLeftEdge;
        vz( this.btmLeftIndecesX, this.btmLeftIndecesY ) = btmLeftEdge;
        vz( this.topRightIndecesX, this.topRightIndecesY ) = topRightEdge;
        vz( this.btmRightIndecesX, this.btmRightIndecesY ) = btmRightEdge;
    end
    
    function [ extrapolateEdge ] = ExtrapolateEdge( this, edge, t, X, Y, idxX, idxY )
        
        c12 = 1-this.c^2;
        bndFun=(c12*X.^2-(Y-t).^2)./(c12*X.^2+(Y-t).^2).^2;
        edgeVec=edge(:); 
        bndFunVec = bndFun(:);
        
        XX = this.X( idxX , idxY );
        YY = this.Y( idxX , idxY );
        mu =  sum(sum(edgeVec.*bndFunVec)) /sum(sum(bndFunVec.^2));
        if( mu < 0 )
            extrapolateEdge = 0*XX;
        else
            extrapolateEdge = mu*(c12*XX.^2-(YY-t).^2)./(c12*XX.^2+(YY-t).^2).^2;
        end
        
        %{
        figure(2)
        mesh( X(1:10,1)', Y(1,1:10), (extrapolateEdge(1:10,1:10) - edge(1:10,1:10))' );
        figure(3)
        mesh( X(1:10,1)', Y(1,1:10), (extrapolateEdge(1:10,1:10))' );
        figure(4)
        mesh( X(1:10,1)', Y(1,1:10), ( edge(1:10,1:10))' );
        %}
    end
  end
  
  methods( Access = private )
      
    function [ this ] = SetIndecesBoundingBox( this, x, y, numberOfBndPoints )
        
        sx = length( x );
        sy = length( y );
        sxSq = numberOfBndPoints * 2;
        sySq = numberOfBndPoints * 2;
        
        %================= left ======================
        this.leftIndecesX = 1:sx;
        this.leftIndecesY = 1:numberOfBndPoints;
        
        %================= right ======================
        this.rightIndecesX = 1:sx;
        this.rightIndecesY = sy-numberOfBndPoints+1:sy;

        %================== top =====================
        this.topIndecesX = 1:numberOfBndPoints;
        this.topIndecesY = 1:sy;

        %================== btm =====================
        this.btmIndecesX = sy-numberOfBndPoints+1:sy;
        this.btmIndecesY = 1:sy;
                
        %================= top left ======================
        this.topLeftIndecesX = 1:numberOfBndPoints;
        this.topLeftIndecesY = 1:numberOfBndPoints;
        
        %================= top right ======================
        this.topRightIndecesX = 1:numberOfBndPoints;
        this.topRightIndecesY = sy-numberOfBndPoints+1:sy;

        %================== btm left =====================
        this.btmLeftIndecesX = sy-numberOfBndPoints+1:sy;
        this.btmLeftIndecesY = 1:numberOfBndPoints;

        %================== btm right =====================
        this.btmRightIndecesX = sy-numberOfBndPoints+1:sy;
        this.btmRightIndecesY = sy-numberOfBndPoints+1:sy;
        
        %==================================================
        %================= top left ======================
        this.topLeftSqIndecesX = numberOfBndPoints+1:sxSq;
        this.topLeftSqIndecesY = numberOfBndPoints+1:sySq;
        
        %================= top right ======================
        this.topRightSqIndecesX = numberOfBndPoints+1:sxSq;
        this.topRightSqIndecesY = sy-sySq+1:sy-numberOfBndPoints;

        %================== btm left =====================
        this.btmLeftSqIndecesX = sy-sxSq+1:sy-numberOfBndPoints;
        this.btmLeftSqIndecesY = numberOfBndPoints+1:sySq;

        %================== btm right =====================
        this.btmRightSqIndecesX = sy-sxSq+1:sy-numberOfBndPoints;
        this.btmRightSqIndecesY = sy-sySq+1:sy-numberOfBndPoints;

    end
    
  end
    
end