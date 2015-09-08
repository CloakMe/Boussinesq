function dy = SIRModel( t, y )  %S, I, R )
 
    dy = zeros(2,1);    % a column vector
    v = y(1);
    u = y(2);   
    
    dy(1) = u;
    dy(2) = -v;
end