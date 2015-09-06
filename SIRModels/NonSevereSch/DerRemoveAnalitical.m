function y = DerRemoveAnalitical( t, S_0, N )
    
    [ a, r, ro ] = GetParams(); 

    q = ( S_0/ro - 1 );
    alphaSqr = ( q^2 + 2 * S_0 * ( N - S_0 )/ro^2 );
    alpha = sqrt( alphaSqr );
    
    p1 = a * alphaSqr * ro^2 /( 2 * S_0 );
    p2 = a*alpha;
    p3 = atanh( q ) / alpha;
    
    y = p1 * sech( p2/2 * t - p3 ).^2;
end