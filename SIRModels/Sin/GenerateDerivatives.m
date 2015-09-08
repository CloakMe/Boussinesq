function ic = GenerateDerivatives( order, IC )
    
    if( order < 0 )
        error('order cannot be negative');
    end
    if( order == 0 || order == 1 )
        ic = IC;
        return
    end
    ic = zeros( order + 1, length( IC ) );
    ic(1,:) = IC;
    for i = 1:order
        ic(i+1,:) = SinModel( 0, ic(i,:) );        
    end

end
