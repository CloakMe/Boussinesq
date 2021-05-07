function [ idx ] = GetIdx( vec, val )
    for k = 1:length(vec)
        if( abs( vec(k) - val ) < 10^(-6) )
            idx=k;return;
        end
    end
    error('could not find the given val incide the vector vec!');
end