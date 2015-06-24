classdef (ConstructOnLoad) BEDomainUtilsP2 < BEDomainUtils
  properties (SetAccess = private, GetAccess = public)
    beta
    c
    mu
    theta
  end
  methods
      
    function this = BEDomainUtilsP2( x, y, order, beta, c, mu, theta )
        this = this@BEDomainUtils( x, y, order );
        this.beta = beta;
        this.c = c;
        this.mu = mu;
        this.theta = theta;
    end
    
    function [ resx ] = TestFun( this, X, Y )
        resx = X.^2 + Y.^2;
    end
    
    function testFunMid = GetTestFunMid( this )
        testFunMid = TestFun( this, this.X, this.Y );
    end
    
    function testFunLeft = GetTestFunLeft( this )
        testFunLeft = TestFun( this, this.leftX, this.leftY );
    end
    
    function testFunRight = GetTestFunRight( this )
        testFunRight = TestFun( this, this.rightX, this.rightY );
    end
    
    function testFunTop = GetTestFunTop( this )
        testFunTop = TestFun( this, this.topX, this.topY );
    end
    
    function testFunBtm = GetTestFunBtm( this )
        testFunBtm = TestFun( this, this.btmX, this.btmY );
    end
       
    function [sdah] = DersBnd( this, X, Y, t, derOrd )
        
        if( nargin == 4 )
            derOrd = this.order;
        end
        sdah = zeros( size( X, 1 ), size( X, 2 ), this.order+1 );
        
        c = this.c;
        bt = this.beta;
        muTheta = this.mu * this.theta;
        %t=t/sqrt(bt);        
        A = (1-c^2)*X.^2 - Y.^2;
        B = (1-c^2)*X.^2 + Y.^2;
        
        sdah(:,:,1) =  muTheta * ( (1-c^2) * X.^2 - (Y - c * t).^2 ) ./ ((1-c^2) * X.^2 + (Y - c * t).^2).^2;
        
        if( derOrd == 0 )
            return;
        end
        
        sdah(:,:,2) =  2 * muTheta * c * (Y - c * t) .* (B + 2 * Y * c * t - c ^ 2 * t ^ 2 + 2 * A) ./...
                       (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 3;
        
        if( derOrd == 1 )
            return;
        end
        
        sdah(:,:,3) =  -2 * muTheta * c ^ 2 * (-10 * A * c ^ 2 * t ^ 2 - 8 * B * c ^ 2 * t ^ 2 -...
                        8 * B .* Y .^ 2 + 16 * Y .^ 2 * c ^ 2 * t ^ 2 - 12 * Y * c ^ 3 * t ^ 3 -...
                        8 * Y .^ 3 * c * t + B .^ 2 + 3 * c ^ 4 * t ^ 4 + 2 * A .* B - 12 * A .* Y .^ 2 +...
                        16 * B .* Y * c * t + 20 * A .* Y * c * t) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 4;
                    
        if( this.order == 2 || derOrd == 2 )
            return;
        end

        sdah(:,:,4) = -24 * muTheta * c ^ 3 * (Y - c * t) .*...
                    (-5 * A * c ^ 2 * t ^ 2 - 5 * B * c ^ 2 * t ^ 2 - 6 * B .* Y .^ 2 +...
                    6 * Y .^ 2 * c ^ 2 * t ^ 2 - 4 * Y * c ^ 3 * t ^ 3 - 4 * Y .^ 3 * c * t +...
                    2 * B .^ 2 + c ^ 4 * t ^ 4 + 3 * A .* B - 8 * A .* Y .^ 2 + 10 * B .* Y * c * t +...
                    10 * A .* Y * c * t) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 5;
                
        if( this.order == 3 || derOrd == 3 )
            return;
        end

        sdah(:,:,5) = 24 * muTheta * c ^ 4 *...
                    (35 * A * c ^ 4 * t ^ 4 + 84 * A .* Y * c * t .* B + 2 * B .^ 3 -...
                    208 * Y .^ 3 .* B * c * t + 264 * Y .^ 2 .* B * c ^ 2 * t ^ 2 +...
                    66 * Y .* c * t .* B .^ 2 - 160 * Y .* c ^ 3 * t ^ 3 .* B - 224 * A .* Y .^ 3 * c * t +...
                    252 * A .* Y .^ 2 * c ^ 2 * t ^ 2 - 140 * A .* Y * c ^ 3 * t ^ 3 -...
                    96 * Y .^ 4 * c ^ 2 * t ^ 2 + 30 * Y * c ^ 5 * t ^ 5 - 80 * Y .^ 2 * c ^ 4 * t ^ 4 +...
                    32 * Y .^ 5 * c * t + 120 * Y .^ 3 * c ^ 3 * t ^ 3 + 64 * Y .^ 4 .* B +...
                    80 * A .* Y .^ 4 - 36 * Y .^ 2 .* B .^ 2 + 3 * A .* B .^ 2 - 5 * c ^ 6 * t ^ 6 +...
                    40 * c ^ 4 * t ^ 4 .* B - 33 * c ^ 2 * t ^ 2 * B .^ 2 -...
                    42 * A .* c ^ 2 * t ^ 2 .* B - 48 * B .* A .* Y .^ 2) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 6;
                
        if( this.order == 4 || derOrd == 4 )
            return;
        end

        sdah(:,:,6) = 240 * muTheta * c ^ 5 * (Y - c * t) .* (28 * A * c ^ 4 * t ^ 4 +...
                    112 * A .* Y * c * t .* B + 9 * B .^ 3 - 224 * Y .^ 3 .* B * c * t +...
                    252 * Y .^ 2 .* B * c ^ 2 * t ^ 2 + 98 * Y .* c * t .* B .^ 2 -...
                    140 * Y .* c ^ 3 * t ^ 3 .* B - 224 * A .* Y .^ 3 * c * t +...
                    224 * A .* Y .^ 2 * c ^ 2 * t ^ 2 - 112 * A .* Y * c ^ 3 * t ^ 3 -...
                    80 * Y .^ 4 * c ^ 2 * t ^ 2 + 18 * Y * c ^ 5 * t ^ 5 -...
                    52 * Y .^ 2 * c ^ 4 * t ^ 4 + 32 * Y .^ 5 * c * t + 88 * Y .^ 3 * c ^ 3 * t ^ 3 +...
                    80 * Y .^ 4 .* B + 96 * A .* Y .^ 4 - 64 * Y .^ 2 .* B .^ 2 + 12 * A .* B .^ 2 -...
                    3 * c ^ 6 * t ^ 6 + 35 * c ^ 4 * t ^ 4 * B - 49 * c ^ 2 * t ^ 2 * B .^ 2 -...
                    56 * A .* c ^ 2 * t ^ 2 .* B - 80 * B .* A .* Y .^ 2) ./ (B - 2 * Y * c * t + c ^ 2 * t ^ 2) .^ 7;
                
        if(this.order == 5 || derOrd == 5 )
            return;
        end

        sdah(:,:,7) = -720 * muTheta * c ^ 6 * (1664 * Y .^ 5 .* B * c * t -...
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
        
        if(this.order == 6 || derOrd == 6 )
            return;
        end
        
        if(this.order > 6)
            error( 'Time derivatives for the boundary function are up to 6th order!  ' );
        end
    end
        
    function testFunMid = GetDersBndMid( this, t, derOrd )
        if( nargin == 2 )
            testFunMid = DersBnd( this, this.X, this.Y, t );
        else
            testFunMid = DersBnd( this, this.X, this.Y, t, derOrd );
        end
    end
    
    function testFunLeft = GetDersBndLeft( this, t, derOrd )
        if( nargin == 2 )
            testFunLeft = DersBnd( this, this.leftX, this.leftY, t );
        else
            testFunLeft = DersBnd( this, this.leftX, this.leftY, t, derOrd );
        end
    end
    
    function testFunRight = GetDersBndRight( this, t, derOrd )
        if( nargin == 2 )
            testFunRight = DersBnd( this, this.rightX, this.rightY, t );
        else
            testFunRight = DersBnd( this, this.rightX, this.rightY, t, derOrd );
        end        
    end
    
    function testFunTop = GetDersBndTop( this, t, derOrd )
        if( nargin == 2 )
            testFunTop = DersBnd( this, this.topX, this.topY, t ); 
        else
            testFunTop = DersBnd( this, this.topX, this.topY, t, derOrd );
        end        
    end
   
    function testFunBtm = GetDersBndBtm( this, t, derOrd )
        if( nargin == 2 )
            testFunBtm = DersBnd( this, this.btmX, this.btmY, t );
        else
            testFunBtm = DersBnd( this, this.btmX, this.btmY, t, derOrd );
        end        
    end
    
  end
  
end