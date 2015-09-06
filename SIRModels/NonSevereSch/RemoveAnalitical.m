function y = RemoveAnalitical( t, S_0, N )
    
    [ a, r, ro ] = GetParams(); 

    q = ( S_0/ro - 1 );
    alphaSqr = ( q^2 + 2 * S_0 * ( N - S_0 )/ro^2 );
    alpha = sqrt( alphaSqr );
    
    %p1 = a * alphaSqr * ro^2 /( 2 * S_0 );
    p2 = a*alpha;
    p3 = atanh( q ) / alpha;
    
    y = q * r^2 / S_0 *( 1 + alpha / q * tanh( p2/2 * t - p3 ) );
end