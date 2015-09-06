function ic = GenerateDerivativesNonSevere( order, IC )
    
    if( order < 0 )
        error('order cannot be negative');
    end
    if( order == 0 )
        ic = IC;
        return
    end
    ic = zeros( order + 1, length( IC ) );
    ic(1,:) = IC;
    for i = 1:order
        ic = GetDerivativeByOrder( ic, i );
    end

end

function ders = GetDerivativeByOrder( ders, order )

    [ a, r, ro ] = GetParams();

    [ S, I, R, dS, dI, dR, ddS, ddI, ddR, dddS, dddI, dddR, ddddS, ddddI, ddddR ] =...
    GetDers( ders );
    
    switch order 
        case 1 
            derSI = S * I;
            derI = I;
        case 2
            derSI = dS * I + S * dI;
            derI = dI;
        case 3
            derSI = ddS * I + 2 * dS * dI + S * ddI;
            derI = ddI;
        case 4
            derSI = dddS * I + 3 * ddS * dI + 3 * dS * ddI + S * dddI;
            derI = dddI;
        case 5
            derSI = ddddS * I + 4 * dddS * dI + 6 * ddS * ddI + 4 * dS * dddI + S * ddddI;
            derI = ddddI;
    end
            
    dnS = - r * derSI;

    dnI = r * derSI - a * derI;

    dnR = a * derI;

    ders(order+1, :) = [ dnS dnI dnR ];
end

function [ S, I, R, dS, dI, dR, ddS, ddI, ddR, dddS, dddI, dddR, ddddS, ddddI, ddddR ] =...
    GetDers( ders )

    sx = size( ders, 1 );
    sy = size( ders, 2 );
    cur_ders = zeros( 5, sy );
    cur_ders(1:sx,:) = ders;
    %=========================
    
    S = cur_ders( 1, 1 );
    I = cur_ders( 1, 2 );
    R = cur_ders( 1, 3 );
        
    dS = cur_ders( 2, 1 );
    dI = cur_ders( 2, 2 );
    dR = cur_ders( 2, 3 );
    
    ddS = cur_ders( 3, 1 );
    ddI = cur_ders( 3, 2 );
    ddR = cur_ders( 3, 3 );
    
    dddS = cur_ders( 4, 1 );
    dddI = cur_ders( 4, 2 );
    dddR = cur_ders( 4, 3 );
    
    ddddS = cur_ders( 5, 1 );
    ddddI = cur_ders( 5, 2 );
    ddddR = cur_ders( 5, 3 );

end
