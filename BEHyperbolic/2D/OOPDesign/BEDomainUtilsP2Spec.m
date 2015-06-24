classdef (ConstructOnLoad) BEDomainUtilsP2Spec < BEDomainUtilsP2
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
        leftDeltaTimeDerivatives
        rightDeltaTimeDerivatives
        topDeltaTimeDerivatives
        btmDeltaTimeDerivatives
  end
  
  methods
      
    function this = BEDomainUtilsP2Spec( x, y, ord, beta, c, mu, theta )
        this = this@BEDomainUtilsP2( x, y, ord, beta, c, mu, theta );
        numberOfBndPoints = ord/2;
        this = SetIndecesBoundingBox( this, x, y, numberOfBndPoints );
    end
    
    function [ this ] = SetDeltaTimeDerivatives( this, v, t ) %( this, X,               Y,            v, t, timeDerOrd )
        this.leftDeltaTimeDerivatives = this.DeltaTimeDersBnd( this.X( this.leftIndecesX, this.leftIndecesY ),...
                                                               this.Y( this.leftIndecesX, this.leftIndecesY ),...
                                                               v( this.leftIndecesX, this.leftIndecesY ),...
                                                               t,...
                                                               this.ord );
        this.rightDeltaTimeDerivatives = this.DeltaTimeDersBnd( this.X( this.rightIndecesX, this.rightIndecesY ),...
                                                                this.Y( this.rightIndecesX, this.rightIndecesY ),...
                                                                v( this.rightIndecesX, this.rightIndecesY ),...
                                                                t,...
                                                                this.ord );
        this.topDeltaTimeDerivatives = this.DeltaTimeDersBnd( this.X( this.topIndecesX, this.topIndecesY ),...
                                                              this.Y( this.topIndecesX, this.topIndecesY ),...
                                                              v( this.topIndecesX, this.topIndecesY ),...
                                                              t,...
                                                              this.ord );
        this.btmDeltaTimeDerivatives = this.DeltaTimeDersBnd( this.X( this.btmIndecesX, this.btmIndecesY ),...
                                                              this.Y( this.btmIndecesX, this.btmIndecesY ),...
                                                              v( this.btmIndecesX, this.btmIndecesY ),...
                                                              t,...
                                                              this.ord );
    end
        
    function [ zeroMatrix ] = YDer( this, M, finiteDiff )
        if( size( finiteDiff, 1 ) > 1 )
            error( 'Finite Difference expected as a row vector but received column vector!' );
        end
        [ sxM, syM ] = size( M );
        sizeFD = length( finiteDiff );

        zeroMatrix = zeros( size( M ) );
        
        numOfBndPoints = this.ord/2;
        for j=numOfBndPoints+1:size(M,2)-numOfBndPoints
            zeroMatrix(:,j) = M(:,j-numOfBndPoints:j+numOfBndPoints)*finiteDiff';
        end

    end
    
    function [ yder ] = YDeriv( this, M, finiteDiff, timeDerOrd  )
        yder = YDer( this, M, finiteDiff );
        numOfBndPoints = this.ord/2;
        yder(:,1:numOfBndPoints) = this.leftDeltaTimeDerivatives(:,:,timeDerOrd+1);
        yder(:,end-numOfBndPoints+1:end) = this.rightDeltaTimeDerivatives(:,:,timeDerOrd+1);
    end
    
    function [ xder ] = XDeriv( this, M, finiteDiff, timeDerOrd  )
        xder = YDer( this, M', finiteDiff )';
        numOfBndPoints = this.ord/2;
        xder(1:numOfBndPoints,:) = this.topDeltaTimeDerivatives(:,:,timeDerOrd+1);
        xder(end-numOfBndPoints+1:end,:) = this.btmDeltaTimeDerivatives(:,:,timeDerOrd+1);
    end
    
    function [ deltaTimeder ] = DeltaTimeDerevative( this, M, finiteDiff, timeDerOrd  )
        deltaTimeder = YDeriv( this, M, finiteDiff, timeDerOrd ) + ...
                       XDeriv( this, M, finiteDiff, timeDerOrd );        
    end
        
    % === time derivatives
    function deltaTimeDerivatives = DeltaTimeDersBnd( this, X, Y, v, t, timeDerOrd )
        
        if( nargin == 4 )
            timeDerOrd = this.ord;
        end
        
        if( max(size(v) ~= size(X)) || max( size(v) ~= size(Y) ) )
            error( 'Function and domain size missmatch!' );
        end
        
        deltaTimeDerivatives = zeros( size( X, 1 ), size( X, 2 ), this.ord+1 );
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
        %maximum = max( max( abs( A ) ) );
        %minimum = min( min( abs( A ) ) );
        %p = ( maximum - minimum )/10000.0;
        %A( abs( A ) < 10 ) = Inf;
        B = omcs*X.^2 + (Y-c*t).^2;
        bndFun = theta * ( A ./ B.^2 );
        
        mu = v ./ bndFun;
        meanMu = mean( mu(:) );
        mu( abs( A ) < 40 ) = meanMu;
        figure(11)
        mesh( X(:,1), Y(1,:), mu' )
        if( sum ( sum ( abs( v ) > 5 ) ) > 0 )
            breakpoint = true;
        end
        if( abs( max( max ( mu ) ) - min( min( mu ) ) ) > .5 )
            breakpoint = true;
        end
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
                
        if( this.ord == 1 || timeDerOrd == 1 )
            return;
        end

        deltaTimeDerivatives(:,:,3) = 120 * theta * mu .* c ^ 2 * (-1 + omcs) .* (-omcs * X .^ 2 + Y .^ 2 - 2 * Y * c * t +...
            c ^ 2 * t ^ 2) .* (c ^ 4 * t ^ 4 - 4 * Y * c ^ 3 * t ^ 3 + 6 * Y .^ 2 * c ^ 2 * t ^ 2 -... 
            14 * omcs * X .^ 2 * c ^ 2 * t ^ 2 + 28 * omcs * X .^ 2 .* Y * c * t - 4 * Y .^ 3 * c * t -...
            14 * omcs * X .^ 2 .* Y .^ 2 + Y .^ 4 + omcs ^ 2 * X .^ 4) ./ (c ^ 2 * t ^ 2 - 2 * Y * c * t +...
            omcs * X .^ 2 + Y .^ 2).^ 6;
                
        if( this.ord == 2 || timeDerOrd == 2 )
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
        
        if(this.ord == 3 || timeDerOrd == 3 )
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
        if(this.ord == 4 || timeDerOrd == 4 )
            return;
        end
        
        if(this.ord > 4)
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
    
  end
  
  methods( Access = private )
      
    function [ this ] = SetIndecesBoundingBox( this, x, y, numberOfBndPoints )
        
        sx = length( x );
        sy = length( y );
        
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

    end
    
  end
    
end