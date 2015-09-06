function [ T, Y, sol ] = SIRModelTaylorNonSevere( t, Tend, IC, order )
    
    t_vec = [ 1 t t^2/2 t^3/6 t^4/24 ];
    t_vec = t_vec( 1:order+1 );
    vecSize = floor( Tend / t ) + 1;
    
    sol = zeros( order+1, size( IC, 2 ), vecSize );
    Y = zeros( vecSize, size( IC, 2 ) );
    
    sol( :, :, 1 ) = IC;
    Y(1,:) = sol( 1, :, 1 );
    iter = 2;
    while( iter <= vecSize )
        solUp = t_vec*sol( :, :, iter-1 ); 
        sol( :, :, iter ) = GenerateDerivativesNonSevere( order, solUp );
        Y(iter,:) = sol( 1, :, iter );
        iter = iter + 1;
    end
    
    T = 0:t:Tend;

end


