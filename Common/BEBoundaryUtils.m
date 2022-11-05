classdef (ConstructOnLoad) BEBoundaryUtils < BEDomainUtils
  properties (SetAccess = private, GetAccess = public)
    alpha
    beta
    c
    mu
    theta
    h
  end
  
  methods (Access = private)
      
	function [sdah] = AsymptoticFunctionSquare( this, X, Y, t, timeDerOrd)
       sdah = zeros( size(X) );
        
        c = this.c;
        beta = this.beta;
        alpha = this.alpha;
        muTheta = this.mu * this.theta;
        omcs = 1 - c^2;
        
        if( timeDerOrd == 0 )  
            
            sdah = 2 * alpha * beta * muTheta ^ 2 * (-c ^ 2 * t ^ 2 + 2 * c * Y * t + omcs * X .^ 2 - Y .^ 2) .^ 2 ./ ...
                (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 4;

        elseif(timeDerOrd == 1 )
            
            sdah = 8 * alpha * beta * muTheta ^ 2 * (-c ^ 2 * t ^ 2 + 2 * c * Y * t + omcs * X .^ 2 - Y .^ 2) * c .* ...
                (-c * t + Y) .* (-c ^ 2 * t ^ 2 + 2 * c * Y * t + 3 * omcs * X .^ 2 - Y .^ 2) ./ ...
                (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 5;
            return;
        elseif(timeDerOrd == 2 )
            
            sdah =  -8 * alpha * beta * muTheta ^ 2 * c ^ 2 * ...
                (-5 * c ^ 6 * t ^ 6 + 30 * c ^ 5 * Y * t ^ 5 + 33 * c ^ 4 * omcs * X .^ 2 * t ^ 4 - ...
                75 * c ^ 4 * Y .^ 2 * t ^ 4 - 132 * c ^ 3 * Y * omcs .* X .^ 2 * t ^ 3 + 100 * c ^ 3 * Y .^ 3 * t ^ 3 - ...
                39 * c ^ 2 * omcs ^ 2 * X .^ 4 * t ^ 2 + 198 * c ^ 2 * omcs * X .^ 2 .* Y .^ 2 * t ^ 2 - ...
                75 * c ^ 2 * Y .^ 4 * t ^ 2 - 132 * c * Y .^ 3 * omcs .* X .^ 2 * t + 78 * c * Y * omcs ^ 2 .* X .^ 4 * t + ...
                30 * c * Y .^ 5 * t - 5 * Y .^ 6 - 39 * Y .^ 2 * omcs ^ 2 .* X .^ 4 + 33 * Y .^ 4 * omcs .* X .^ 2 + ...
                3 * omcs ^ 3 .* X .^ 6) ./ (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 6;
            return;
        elseif(timeDerOrd == 3 )
            
            sdah = -48 * alpha * beta * muTheta ^ 2 * c ^ 3 * (-c * t + Y) .* ...
                (-5 * c ^ 6 * t ^ 6 + 30 * c ^ 5 * Y * t ^ 5 + 49 * c ^ 4 * omcs * X .^ 2 * t ^ 4 - ...
                75 * c ^ 4 * Y .^ 2 * t ^ 4 - 196 * c ^ 3 * Y * omcs .* X .^ 2 * t ^ 3 + 100 * c ^ 3 * Y .^ 3 * t ^ 3 - ...
                87 * c ^ 2 * omcs ^ 2 * X .^ 4 * t ^ 2 + 294 * c ^ 2 * omcs * X .^ 2 .* Y .^ 2 * t ^ 2 - ...
                75 * c ^ 2 * Y .^ 4 * t ^ 2 - 196 * c * Y .^ 3 * omcs .* X .^ 2 * t + 174 * c * Y * omcs ^ 2 .* X .^ 4 * t + ...
                30 * c * Y .^ 5 * t - 5 * Y .^ 6 - 87 * Y .^ 2 * omcs ^ 2 .* X .^ 4 + ...
                49 * Y .^ 4 * omcs .* X .^ 2 + 19 * omcs ^ 3 * X .^ 6) ./ ...
                (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 7;
            return;
        elseif(timeDerOrd == 4 )
            
            sdah = 48 * alpha * beta * muTheta ^ 2 * c ^ 4 * (35 * c ^ 8 * t ^ 8 - 280 * c ^ 7 * Y * t ^ 7 - ...
                476 * c ^ 6 * omcs * X .^ 2 * t ^ 6 + 980 * c ^ 6 * Y .^ 2 * t ^ 6 - 1960 * c ^ 5 * Y .^ 3 * t ^ 5 + ...
                2856 * c ^ 5 * Y * omcs .* X .^ 2 * t ^ 5 + 1202 * c ^ 4 * omcs ^ 2 * X .^ 4 * t ^ 4 + ...
                2450 * c ^ 4 * Y .^ 4 * t ^ 4 - 7140 * c ^ 4 * omcs * X .^ 2 .* Y .^ 2 * t ^ 4 + ...
                9520 * c ^ 3 * Y .^ 3 * omcs .* X .^ 2 * t ^ 3 - 1960 * c ^ 3 * Y .^ 5 * t ^ 3 - ...
                4808 * c ^ 3 * Y * omcs ^ 2 .* X .^ 4 * t ^ 3 + 980 * c ^ 2 * Y .^ 6 * t ^ 2 - ...
                508 * c ^ 2 * omcs ^ 3 * X .^ 6 * t ^ 2 - 7140 * c ^ 2 * Y .^ 4 * omcs .* X .^ 2 * t ^ 2 + ...
                7212 * c ^ 2 * Y .^ 2 * omcs ^ 2 .* X .^ 4 * t ^ 2 - 280 * c * Y .^ 7 * t + 2856 * c * Y .^ 5 * omcs .* X .^ 2 * t + ...
                1016 * c * Y * omcs ^ 3 .* X .^ 6 * t - 4808 * c * Y .^ 3 * omcs ^ 2 .* X .^ 4 * t - ...
                508 * omcs ^ 3 .* X .^ 6 .* Y .^ 2 - 476 * omcs .* X .^ 2 .* Y .^ 6 + 35 * Y .^ 8 + ...
                1202 * omcs ^ 2 * X .^ 4 .* Y .^ 4 + 19 * omcs ^ 4 .* X .^ 8) ./ ...
                (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 8;
            return;
        end
    end
    
	function [sdah] = FPSOperatorAsymptoticFunction( this, X, Y, t, timeDerOrd)
       sdah = zeros( size(X) );
        
        c = this.c;
        beta = this.beta;
        muTheta = this.mu * this.theta;
        omcs = 1 - c^2;
        
        if( timeDerOrd == 0 )  
            
            sdah = -6 * muTheta * c ^ 2 * (c ^ 4 * t ^ 4 - 4 * c ^ 3 * Y * t ^ 3 - 6 * c ^ 2 * omcs * X .^ 2 * t ^ 2 +...
                6 * c ^ 2 * Y .^ 2 * t ^ 2 + 12 * c * Y * omcs .* X .^ 2 * t - 4 * c * Y .^ 3 * t -...
                6 * omcs * X .^ 2 .* Y .^ 2 + omcs ^ 2 * X .^ 4 + Y .^ 4) ./ ...
                (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 4;

            return;
        elseif(timeDerOrd == 1 )
            
            sdah = -24 * muTheta * c ^ 3 * (-c * t + Y) .* (c ^ 4 * t ^ 4 - 4 * c ^ 3 * Y * t ^ 3 - 10 * c ^ 2 * omcs * X .^ 2 * t ^ 2 +...
                6 * c ^ 2 * Y .^ 2 * t ^ 2 + 20 * c * Y * omcs .* X .^ 2 * t - 4 * c * Y .^ 3 * t - 10 * omcs * X .^ 2 .* Y .^ 2 +...
                5 * omcs ^ 2 * X .^ 4 + Y .^ 4) ./ (c ^ 2 * t ^ 2 - 2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 5;

            return;
        elseif(timeDerOrd == 2 )
            
            sdah =  120 * muTheta * c ^ 4 * (omcs * X .^ 2 - Y .^ 2 + 2 * c * Y * t - c ^ 2 * t ^ 2) .* (omcs ^ 2 * X .^ 4 -...
                14 * omcs * X .^ 2 .* Y .^ 2 - 14 * c ^ 2 * omcs * X .^ 2 * t ^ 2 + 28 * c * Y * omcs .* X .^ 2 * t + Y .^ 4 +...
                6 * c ^ 2 * Y .^ 2 * t ^ 2 + c ^ 4 * t ^ 4 - 4 * c ^ 3 * Y * t ^ 3 - 4 * c * Y .^ 3 * t) ./ (c ^ 2 * t ^ 2 -...
                2 * c * Y * t + omcs * X .^ 2 + Y .^ 2) .^ 6;

            return;
        elseif(timeDerOrd == 3 )
            
            sdah = 720 * muTheta * c ^ 5 * (-c * t + Y) .* (-c ^ 6 * t ^ 6 + 6 * c ^ 5 * Y * t ^ 5 +...
                21 * omcs * X .^ 2 * c ^ 4 * t ^ 4 - 15 * Y .^ 2 * c ^ 4 * t ^ 4 + 20 * Y .^ 3 * c ^ 3 * t ^ 3 -...
                84 * omcs * X .^ 2 * c ^ 3 .* Y * t ^ 3 + 126 * omcs * X .^ 2 .* Y .^ 2 * c ^ 2 * t ^ 2 -...
                15 * Y .^ 4 * c ^ 2 * t ^ 2 - 35 * omcs ^ 2 * X .^ 4 * c ^ 2 * t ^ 2 - 84 * omcs * X .^ 2 .* Y .^ 3 * c * t +...
                6 * Y .^ 5 * c * t + 70 * omcs ^ 2 * X .^ 4 * c .* Y * t + 7 * omcs ^ 3 * X .^ 6 - Y .^ 6 - ...
                35 * omcs ^ 2 * X .^ 4 .* Y .^ 2 + 21 * omcs * X .^ 2 .* Y .^ 4) ./ (c ^ 2 * t ^ 2 - 2 * c * Y * t +...
                omcs * X .^ 2 + Y .^ 2) .^ 7;

            return;
        elseif(timeDerOrd == 4 )
            
            sdah = -5040 * muTheta * c ^ 6 * (c ^ 8 * t ^ 8 - 8 * c ^ 7 * Y * t ^ 7 - 28 * c ^ 6 * omcs * X .^ 2 * t ^ 6 +...
                28 * c ^ 6 * Y .^ 2 * t ^ 6 + 168 * c ^ 5 * Y * omcs .* X .^ 2 * t ^ 5 - 56 * c ^ 5 * Y .^ 3 * t ^ 5 +...
                70 * c ^ 4 * Y .^ 4 * t ^ 4 - 420 * c ^ 4 * omcs * X .^ 2 .* Y .^ 2 * t ^ 4 + 70 * c ^ 4 * omcs ^ 2 .* X .^ 4 * t ^ 4 -....
                56 * c ^ 3 * Y .^ 5 * t ^ 3 + 560 * c ^ 3 * Y .^ 3 * omcs .* X .^ 2 * t ^ 3 - 280 * c ^ 3 * Y * omcs ^ 2 .* X .^ 4 * t ^ 3 +...
                420 * c ^ 2 * Y .^ 2 * omcs ^ 2 .* X .^ 4 * t ^ 2 - 28 * c ^ 2 * omcs ^ 3 * X .^ 6 * t ^ 2 +...
                28 * c ^ 2 * Y .^ 6 * t ^ 2 - 420 * c ^ 2 * Y .^ 4 * omcs .* X .^ 2 * t ^ 2 -...
                280 * c * Y .^ 3 * omcs ^ 2 .* X .^ 4 * t + 168 * c * Y .^ 5 * omcs .* X .^ 2 * t - 8 * c * Y .^ 7 * t +...
                56 * c * Y * omcs ^ 3 .* X .^ 6 * t + Y .^ 8 + omcs ^ 4 * X .^ 8 + 70 * Y .^ 4 * omcs ^ 2 .* X .^ 4 -...
                28 * Y .^ 6 * omcs .* X .^ 2 - 28 * Y .^ 2 * omcs ^ 3 .* X .^ 6) ./ (c ^ 2 * t ^ 2 - 2 * c * Y * t +...
                omcs * X .^ 2 + Y .^ 2) .^ 8;

            return;
        end
    end
    
    function [sdah] = AsymptoticFunction( this, X, Y, t, timeDerOrd )
        
        sdah = zeros( size(X) );
        
        c = this.c;
        bt = this.beta;
        muTheta = this.mu * this.theta;
        %t=t/sqrt(bt);        
        A = (1-c^2)*X.^2 - Y.^2;
        B = (1-c^2)*X.^2 + Y.^2;
                
        if( timeDerOrd == 0 )
            sdah =  muTheta * ( (1-c^2) * X.^2 - (Y - c * t).^2 ) ./ ((1-c^2) * X.^2 + (Y - c * t).^2).^2;
            return;
        elseif(timeDerOrd == 1 )
        
            sdah =  2 * muTheta * c * (Y - c * t) .* (B + 2 * Y * c * t - c ^ 2 * t ^ 2 + 2 * A) ./...
                       (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 3;
            return;
        elseif(timeDerOrd == 2 )
        
            sdah =  -2 * muTheta * c ^ 2 * (-10 * A * c ^ 2 * t ^ 2 - 8 * B * c ^ 2 * t ^ 2 -...
                        8 * B .* Y .^ 2 + 16 * Y .^ 2 * c ^ 2 * t ^ 2 - 12 * Y * c ^ 3 * t ^ 3 -...
                        8 * Y .^ 3 * c * t + B .^ 2 + 3 * c ^ 4 * t ^ 4 + 2 * A .* B - 12 * A .* Y .^ 2 +...
                        16 * B .* Y * c * t + 20 * A .* Y * c * t) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 4;
                    
            return;
        elseif(timeDerOrd == 3 )

        	sdah = -24 * muTheta * c ^ 3 * (Y - c * t) .*...
                    (-5 * A * c ^ 2 * t ^ 2 - 5 * B * c ^ 2 * t ^ 2 - 6 * B .* Y .^ 2 +...
                    6 * Y .^ 2 * c ^ 2 * t ^ 2 - 4 * Y * c ^ 3 * t ^ 3 - 4 * Y .^ 3 * c * t +...
                    2 * B .^ 2 + c ^ 4 * t ^ 4 + 3 * A .* B - 8 * A .* Y .^ 2 + 10 * B .* Y * c * t +...
                    10 * A .* Y * c * t) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 5;
                
            return;
        elseif(timeDerOrd == 4 )

        	sdah = 24 * muTheta * c ^ 4 *...
                    (35 * A * c ^ 4 * t ^ 4 + 84 * A .* Y * c * t .* B + 2 * B .^ 3 -...
                    208 * Y .^ 3 .* B * c * t + 264 * Y .^ 2 .* B * c ^ 2 * t ^ 2 +...
                    66 * Y .* c * t .* B .^ 2 - 160 * Y .* c ^ 3 * t ^ 3 .* B - 224 * A .* Y .^ 3 * c * t +...
                    252 * A .* Y .^ 2 * c ^ 2 * t ^ 2 - 140 * A .* Y * c ^ 3 * t ^ 3 -...
                    96 * Y .^ 4 * c ^ 2 * t ^ 2 + 30 * Y * c ^ 5 * t ^ 5 - 80 * Y .^ 2 * c ^ 4 * t ^ 4 +...
                    32 * Y .^ 5 * c * t + 120 * Y .^ 3 * c ^ 3 * t ^ 3 + 64 * Y .^ 4 .* B +...
                    80 * A .* Y .^ 4 - 36 * Y .^ 2 .* B .^ 2 + 3 * A .* B .^ 2 - 5 * c ^ 6 * t ^ 6 +...
                    40 * c ^ 4 * t ^ 4 .* B - 33 * c ^ 2 * t ^ 2 * B .^ 2 -...
                    42 * A .* c ^ 2 * t ^ 2 .* B - 48 * B .* A .* Y .^ 2) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 6;
                
            return;
        elseif(timeDerOrd == 5 )

            sdah = 240 * muTheta * c ^ 5 * (Y - c * t) .* (28 * A * c ^ 4 * t ^ 4 +...
                    112 * A .* Y * c * t .* B + 9 * B .^ 3 - 224 * Y .^ 3 .* B * c * t +...
                    252 * Y .^ 2 .* B * c ^ 2 * t ^ 2 + 98 * Y .* c * t .* B .^ 2 -...
                    140 * Y .* c ^ 3 * t ^ 3 .* B - 224 * A .* Y .^ 3 * c * t +...
                    224 * A .* Y .^ 2 * c ^ 2 * t ^ 2 - 112 * A .* Y * c ^ 3 * t ^ 3 -...
                    80 * Y .^ 4 * c ^ 2 * t ^ 2 + 18 * Y * c ^ 5 * t ^ 5 -...
                    52 * Y .^ 2 * c ^ 4 * t ^ 4 + 32 * Y .^ 5 * c * t + 88 * Y .^ 3 * c ^ 3 * t ^ 3 +...
                    80 * Y .^ 4 .* B + 96 * A .* Y .^ 4 - 64 * Y .^ 2 .* B .^ 2 + 12 * A .* B .^ 2 -...
                    3 * c ^ 6 * t ^ 6 + 35 * c ^ 4 * t ^ 4 * B - 49 * c ^ 2 * t ^ 2 * B .^ 2 -...
                    56 * A .* c ^ 2 * t ^ 2 .* B - 80 * B .* A .* Y .^ 2) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 7;
                
            return;
        elseif(timeDerOrd == 6 )

            sdah = -720 * muTheta * c ^ 6 * (1664 * Y .^ 5 .* B * c * t -...
                3072 * Y .^ 4 .* B * c ^ 2 * t ^ 2 - 1264 * Y .^ 3 * c * t .* B .^ 2 +...
                3136 * Y .^ 3 * c ^ 3 * t ^ 3 .* B + 1584 * Y .^ 2 * c ^ 2 * t ^ 2 .* B .^ 2 -...
                1904 * Y .^ 2 * c ^ 4 * t ^ 4 .* B - 952 * Y .* c ^ 3 * t ^ 3 .* B .^ 2 +...
                672 * Y .* c ^ 5 * t ^ 5 .* B + 176 * Y .* c * t .* B .^ 3 + 1728 * A .* Y .^ 5 * c * t -...
                2880 * A .* Y .^ 4 * c ^ 2 * t ^ 2 + 2688 * A .* Y .^ 3 * c ^ 3 * t ^ 3 -...
                1512 * A .* Y .^ 2 * c ^ 4 * t ^ 4 + 504 * A .* Y * c ^ 5 * t ^ 5 +...
                252 * A .* c ^ 4 * t ^ 4 .* B - 108 * A .* c ^ 2 * t ^ 2 .* B .^ 2 +...
                3 * B .^ 4 - 1440 * A .* Y .^ 3 * c * t .* B + 1728 * A .* Y .^ 2 * c ^ 2 * t ^ 2 .* B -...
                1008 * A .* Y * c ^ 3 * t ^ 3 .* B + 216 * A .* Y * c * t .* B .^ 2 +...
                4 * B .^ 3 .* A + 400 * Y .^ 4 .* B .^ 2 + 7 * c ^ 8 * t ^ 8 -...
                96 * Y .^ 2 .* B .^ 3 - 384 * Y .^ 6 .* B - 448 * A .* Y .^ 6 -...
                84 * A * c ^ 6 * t ^ 6 + 238 * c ^ 4 * t ^ 4 * B .^ 2 - 112 * c ^ 6 * t ^ 6 * B +...
                896 * Y .^ 4 * c ^ 4 * t ^ 4 - 560 * Y .^ 3 * c ^ 5 * t ^ 5 -...
                120 * B .^ 2 .* A .* Y .^ 2 - 88 * c ^ 2 * t ^ 2 * B .^ 3 -...
                896 * Y .^ 5 * c ^ 3 * t ^ 3 - 56 * Y * c ^ 7 * t ^ 7 + 480 * A .* Y .^ 4 .* B -...
                128 * Y .^ 7 * c * t + 512 * Y .^ 6 * c ^ 2 * t ^ 2 + 224 * Y .^ 2 * c ^ 6 * t ^ 6) ./...
                (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 8;
        
            return;
        else
            error( 'Time derivatives for the boundary function are up to 6th order!  ' );
        end
     end
  end
  
  methods (Access = public)
      
    function this = BEBoundaryUtils( x, y, order, alpha, beta, c, mu, theta, h )
        this = this@BEDomainUtils( x, y, order);
        this.alpha = alpha;
        this.beta = beta;
        this.c = c;
        this.mu = mu;
        this.theta = theta;
        this.h = h;
    end

    function [ result ] = AssymptoticFuncSquareOutside( this, finiteDiff, t, timeDerOrder ) 

        bndPntsCount = this.order/2;
        vdah = zeros(length(this.x), length(this.y));
        
        xNet = this.X_xAugDomain(1:bndPntsCount,:);
        yNet = this.Y_xAugDomain(1:bndPntsCount,:);
        asymptoticFunctionXStart = this.AsymptoticFunctionSquare( xNet, yNet, t, timeDerOrder );
        xNet = this.X_xAugDomain(end - bndPntsCount + 1:end,:);
        yNet = this.Y_xAugDomain(end - bndPntsCount + 1:end,:);
        asymptoticFunctionXEnd = this.AsymptoticFunctionSquare( xNet, yNet, t, timeDerOrder );

        clear('xNet'); clear('yNet');

        xNet = this.X_yAugDomain(:,1:bndPntsCount);
        yNet = this.Y_yAugDomain(:,1:bndPntsCount);
        asymptoticFunctionYStart = AsymptoticFunctionSquare( this, xNet, yNet, t, timeDerOrder );
        xNet = this.X_yAugDomain(:,end - bndPntsCount + 1:end);
        yNet = this.Y_yAugDomain(:,end - bndPntsCount + 1:end);
        asymptoticFunctionYEnd = AsymptoticFunctionSquare( this, xNet, yNet, t, timeDerOrder );
               
        result = (...
             this.YDerivative( vdah, asymptoticFunctionYStart, asymptoticFunctionYEnd, finiteDiff ) +...
             this.XDerivative( vdah, asymptoticFunctionXStart, asymptoticFunctionXEnd, finiteDiff ));      %/this.h^2           
         
        %figure(1)
        %yy = (this.y(1)-3*this.h):this.h:(this.y(1)+3*this.h);
        %mesh(this.x', yy' , result(:,1:7)');
        %xlabel('x');            ylabel('y');
        %yoyo = 6;
    end
        
    function [ result ] = FPSOperatorOutside( this, finiteDiff, t, timeDerOrder ) 

        bndPntsCount = this.order/2;
        vdah = zeros(length(this.x), length(this.y));
        
        xNet = this.X_xAugDomain(1:bndPntsCount,:);
        yNet = this.Y_xAugDomain(1:bndPntsCount,:);
        asymptoticFunctionXStart = this.FPSOperatorAsymptoticFunction( xNet, yNet, t, timeDerOrder );
        xNet = this.X_xAugDomain(end - bndPntsCount + 1:end,:);
        yNet = this.Y_xAugDomain(end - bndPntsCount + 1:end,:);
        asymptoticFunctionXEnd = this.FPSOperatorAsymptoticFunction( xNet, yNet, t, timeDerOrder );

        clear('xNet'); clear('yNet');

        xNet = this.X_yAugDomain(:,1:bndPntsCount);
        yNet = this.Y_yAugDomain(:,1:bndPntsCount);
        asymptoticFunctionYStart = FPSOperatorAsymptoticFunction( this, xNet, yNet, t, timeDerOrder );
        xNet = this.X_yAugDomain(:,end - bndPntsCount + 1:end);
        yNet = this.Y_yAugDomain(:,end - bndPntsCount + 1:end);
        asymptoticFunctionYEnd = FPSOperatorAsymptoticFunction( this, xNet, yNet, t, timeDerOrder );
               
        result = (...
             this.YDerivative( vdah, asymptoticFunctionYStart, asymptoticFunctionYEnd, finiteDiff ) +...
             this.XDerivative( vdah, asymptoticFunctionXStart, asymptoticFunctionXEnd, finiteDiff ));      %/this.h^2           
         
        %figure(1)
        %yy = (this.y(1)-3*this.h):this.h:(this.y(1)+3*this.h);
        %mesh(this.x', yy' , result(:,1:7)');
        %xlabel('x');            ylabel('y');
        %yoyo = 6;
    end
    
    function [ result ] = AssymptoticFuncOutside( this, coeff, finiteDiff, t, timeDerOrder ) 

        bndPntsCount = this.order/2;
        vdah = zeros(length(this.x), length(this.y));
        dnU_dtnOnBoundary = vdah;
        
        if(timeDerOrder > 0)
            fdStartPos = this.order/2;
            timeDerivative = BEUtilities.GetFinDiffCoeff( -fdStartPos:fdStartPos, timeDerOrder)';
            xNet = this.X_yAugDomain(:,1:bndPntsCount);
            yNet = this.Y_yAugDomain(:,1:bndPntsCount);
            asymptoticFunctionYStart = this.AsymptoticFunction( xNet, yNet, t, timeDerOrder );      
            xNet = this.X_yAugDomain(:,end - bndPntsCount + 1:end);
            yNet = this.Y_yAugDomain(:,end - bndPntsCount + 1:end);
            asymptoticFunctionYEnd = this.AsymptoticFunction( xNet, yNet, t, timeDerOrder );
            %timeDerivative could be first, second, third or fourth
            %derivative! finiteDiff is spatial second derivative.
            h = this.y(2) - this.y(1);
            dnU_dtnOnBoundary = (-this.c)^timeDerOrder * this.YDerivative( vdah, asymptoticFunctionYStart, asymptoticFunctionYEnd, timeDerivative )/h^timeDerOrder;  
            clear('xNet'); clear('yNet');
        end
        
        xNet = this.X_xAugDomain(1:bndPntsCount,:);
        yNet = this.Y_xAugDomain(1:bndPntsCount,:);
        asymptoticFunctionXStart = this.AsymptoticFunction( xNet, yNet, t, timeDerOrder );
        xNet = this.X_xAugDomain(end - bndPntsCount + 1:end,:);
        yNet = this.Y_xAugDomain(end - bndPntsCount + 1:end,:);
        asymptoticFunctionXEnd = this.AsymptoticFunction( xNet, yNet, t, timeDerOrder );

        clear('xNet'); clear('yNet');

        xNet = this.X_yAugDomain(:,1:bndPntsCount);
        yNet = this.Y_yAugDomain(:,1:bndPntsCount);
        asymptoticFunctionYStart = AsymptoticFunction( this, xNet, yNet, t, timeDerOrder );
        xNet = this.X_yAugDomain(:,end - bndPntsCount + 1:end);
        yNet = this.Y_yAugDomain(:,end - bndPntsCount + 1:end);
        asymptoticFunctionYEnd = AsymptoticFunction( this, xNet, yNet, t, timeDerOrder );
               
        %timeDerivative inside d^nU/dt^n could be first, second, third or fourth
        %derivative! finiteDiff is spatial second derivative.
        result = coeff*(...
             dnU_dtnOnBoundary + ...
             this.YDerivative( vdah, asymptoticFunctionYStart, asymptoticFunctionYEnd, finiteDiff ) +...
             this.XDerivative( vdah, asymptoticFunctionXStart, asymptoticFunctionXEnd, finiteDiff ));                 

        %if(nargin == 7)
        %    deltaV = this.beta*this.DeltaH(finiteDiff,v)/this.h^2;
        %    newV = [deltaAsymptoticFunctionYStart deltaV deltaAsymptoticFunctionYEnd];
        %    figure(1)
        %    yy = this.y(1)-3*this.h:h:this.y(1)+3*this.h;
        %    plot(this.x, yy , newV(1:7,:)');
        %end
        
    end
    
    function [ dudt_t0 ] = dudtInside_t0( this, finiteDiff, u_t0 ) 

        %dersBndLeft = domUtils.GetAsymptoticFunctionLeft( 0 );
            %dersBndRight = domUtils.GetDersBndRight( 0 );
            bndPntsCount = this.order/2;
            xNet = this.X_yAugDomain(:,1:bndPntsCount);
            yNet = this.Y_yAugDomain(:,1:bndPntsCount);
            asymptoticFunctionYStart = this.AsymptoticFunction( xNet, yNet, 0, 0 );      
            xNet = this.X_yAugDomain(:,end - bndPntsCount + 1:end);
            yNet = this.Y_yAugDomain(:,end - bndPntsCount + 1:end);
            asymptoticFunctionYEnd = this.AsymptoticFunction( xNet, yNet, 0, 0 );

            dudt_t0 =  -this.c * this.YDerivative( u_t0, asymptoticFunctionYStart, asymptoticFunctionYEnd, finiteDiff )/this.h;
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
    
    
    function testFunMid = GetAsymptoticFunctionMid( this, t, derOrd )
        if( nargin == 2 )
            testFunMid = AsymptoticFunction( this, this.X, this.Y, t );
        else
            testFunMid = AsymptoticFunction( this, this.X, this.Y, t, derOrd );
        end
    end
    
    function testFunLeft = GetAsymptoticFunctionLeft( this, t, derOrd )
        if( nargin == 2 )
            testFunLeft = AsymptoticFunction( this, this.leftX, this.leftY, t );
        else
            testFunLeft = AsymptoticFunction( this, this.leftX, this.leftY, t, derOrd );
        end
    end
    
    function testFunRight = GetAsymptoticFunctionRight( this, t, derOrd )
        if( nargin == 2 )
            testFunRight = DeltaAsymptoticFunction( this, this.rightX, this.rightY, t );
        else
            testFunRight = DeltaAsymptoticFunction( this, this.rightX, this.rightY, t, derOrd );
        end        
    end
    
    function testFunTop = GetAsymptoticFunctionTop( this, t, derOrd )
        if( nargin == 2 )
            testFunTop = AsymptoticFunction( this, this.topX, this.topY, t ); 
        else
            testFunTop = AsymptoticFunction( this, this.topX, this.topY, t, derOrd );
        end        
    end
    
    function testFunBtm = GetAsymptoticFunctionBtm( this, t, derOrd )
        if( nargin == 2 )
            testFunBtm = AsymptoticFunction( this, this.btmX, this.btmY, t );
        else
            testFunBtm = AsymptoticFunction( this, this.btmX, this.btmY, t, derOrd );
        end        
    end
    
  end
  
end