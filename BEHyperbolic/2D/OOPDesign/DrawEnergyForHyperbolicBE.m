function DrawEnergyForHyperbolicBE( engine, tt )

    en = zeros( size(tt) ); 
    for j = 1:size(tt,2)   
        size( engine.x );
        lenX = length( engine.x );
        lenXHalf = (lenX - 1 )/2;
        lenY = length( engine.y );
        lenYHalf = (lenY - 1 )/2;
        load(['SOL\vu_' num2str(tt(j)) '.mat']);
        size( vu );
        [ crntSx, crntSy ] = size( vu );
        ij = ( crntSx + 1 )/2;
        lo = ( crntSy + 1 )/2;
        vvu = vu( ij - lenXHalf :ij + lenXHalf,...
                  lo - lenYHalf :lo + lenYHalf );
              
        load(['SOL\vz_' num2str(tt(j)) '.mat']);
        vvz = vz( ij - lenXHalf :ij + lenXHalf,...
                  lo - lenYHalf :lo + lenYHalf );
        en( j ) = engine.GetEnergy( vvz, vvu, tt(j) );

        clear('v2u');clear('vvu'); clear('vvz');         
    end
    
    figure( 10 )
    title('Fixed Energy functional');
    xlabel('time "t"');  ylabel('EN');
    plot( tt, en );
    
end