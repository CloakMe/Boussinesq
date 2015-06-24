function [ v ] = BEHeatMassTransferEquation( x, y, n, k, C1, C2 )

    if( nargin < 6 )
        C2 = 0;
    end    
    if( nargin < 5 )
        C1 = 0;
    end    
    if( nargin < 4 )
        k = 1;
    end    
    if( nargin < 3 )
        n = 2;
    end
    
    s = ( ( k * ( 1 - n )^2 ) ^ ( 1/( 1 - n ) ) )/ 4;
    v = s*( ( x + C1 ).^2 + ( y + C2 ).^2 ) .^ ( 1/( 1 - n ) );
        
end
