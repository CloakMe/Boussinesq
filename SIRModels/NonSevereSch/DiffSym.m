clear;
syms t alpha S_0 N ro a r q p1 p2  p3  alphaS 

alpha = sqrt( q^2 + 2 * S_0 * ( N - S_0 )/ro^2 )  ;
alphaS = q^2 + 2 * S_0 * ( N - S_0 )/ro^2  ;
ro = ( a / r )  ;
q = ( S_0/ro - 1 )  ;

p1 = a * alphaS * ro^2/( 2*S_0 );
p2 = a * alpha ;
p3 = atanh( q ) / alpha ;

ires = int( p1 * sech( p2/2 * t - p3 )^2 , t )

(2*(- 2*S_0^2*a^2 + 2*N*S_0*a^2 + a^2*q^2*ro^2)) /...
    (S_0*r^2*ro^2*((- 2*S_0^2 + 2*N*S_0 + q^2*ro^2)/ro^2)^(1/2)*...
    (exp((2*atanh((S_0*r)/a - 1))/(q^2 - (2*S_0^2)/ro^2 + (2*N*S_0)/ro^2)^(1/2) -...
    a*t*(q^2 - (2*S_0^2)/ro^2 + (2*N*S_0)/ro^2)^(1/2)) + 1) )

simp = simplify ( (2*a^2*(- 2*S_0^2 + 2*N*S_0 + q^2*ro^2))/...
    (S_0*r^2*ro^2*...
    ((- 2*S_0^2 + 2*N*S_0 + q^2*ro^2)/ro^2)^(1/2)) )

(- 2*S_0^2*r^2 - 4*S_0*a*r + 4*N*S_0*r^2 + 2*a^2)/(S_0*r^2*((- S_0^2*r^2 - 2*S_0*a*r + 2*N*S_0*r^2 + a^2)/a^2)^(1/2))

return;
clear;
syms t alpha S_0 N ro a r q

alpha = sqrt( q^2 + 2 * S_0 * ( N - S_0 )/ro^2 )  ;
ro = ( a / r )  ;
q = ( S_0/ro - 1 )  ;

res = ...
    diff( r^2/S_0 * ( q + alpha * tanh( a*alpha/2 * t - atanh( q ) / alpha ) ), t )

( a*r^2*alpha^2/(2*S_0)*(tanh(atanh( (S_0*r)/a - 1)/(q^2 + (2*S_0*(N - S_0))/ro^2)^(1/2) -...
    (a*t*(q^2 + (2*S_0*(N - S_0))/ro^2)^(1/2))/2)^2 - 1)    )

