function dy = SIRModelNonSevere( t, y )  %S, I, R )
    
    [ a, r, ro ] = GetParams(); 
    
    dy = zeros(3,1);    % a column vector
    S = y(1);
    I = y(2);
    R = y(3);   
    
    dS = - r * S * I;
    dy(1) = dS;

    dI = r * S * I - a * I;
    dy(2) = dI;

    dR = a * I;
    dy(3) = dR;

end