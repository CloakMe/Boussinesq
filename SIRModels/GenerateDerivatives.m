function ic = GenerateIC( order, IC )
    
    if( order < 0 )
        error('order cannot be negative');
    end
    if( order == 0 )
        ic = IC;
    end
    ic = zeros( order, length( IC ) );
    ic(1,:) = IC;
    for i = 2:order+1
        ic(i,:) = SIRModel( 0, ic(i-1,:) );
    end

end