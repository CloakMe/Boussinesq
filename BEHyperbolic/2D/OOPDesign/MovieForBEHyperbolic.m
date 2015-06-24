function MovieForBEHyperbolic( compBoxToShow, viewTypeX, viewTypeY, tt, x, y )
    mag = 120;
    %mag2=mag;
    idxCompBox = compBoxToShow;
    figure(6)
    clear('F')
    camva('manual');
    %[ ijj ] = GetIdx( x, 0 )
    %[ loo ] = GetIdx( y, 0 );
    [ idxCompBox.x_st ] = GetIdx( x, compBoxToShow.x_st );
    [ idxCompBox.x_end ] = GetIdx( x, compBoxToShow.x_end );
    [ idxCompBox.y_st ] = GetIdx( y, compBoxToShow.y_st );
    [ idxCompBox.y_end ] = GetIdx( y, compBoxToShow.y_end );

    %tt = 1:0.1:20;    
    for j = 1:size(tt,2)   

        xxx=x( idxCompBox.x_st:idxCompBox.x_end );
        lenX = length( xxx );
        lenXHalf = (lenX - 1 )/2;
        yyy=y( idxCompBox.y_st:idxCompBox.y_end ); 
        lenY = length( yyy );
        lenYHalf = (lenY - 1 )/2;
        load(['SOL\vu_' num2str(tt(j)) '.mat']);
        [ crntSx, crntSy ] = size( vu );
        ij = ( crntSx + 1 )/2;
        lo = ( crntSy + 1 )/2;
        %vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
        vvu = vu( ij - lenXHalf :ij + lenXHalf,...
                  lo - lenYHalf :lo + lenYHalf );
        mesh(xxx,yyy,vvu');
        %mesh(x,y,vu');
        title(['Solution at time: ',num2str(tt(j))]);
        xlabel('x'); ylabel('y');
        
        %view(1,90);
        %axis([xxx(1) xxx(end) yyy(1) yyy(end) -.7, 2.5]);
        colorbar;
        %caxis([-.7, 1.7]);
        view( viewTypeX, viewTypeY );
        
        F(:,j) = getframe;
        j=j+1;
        clear('v2u');clear('vu');        
    end
end
function [ idx ] = GetIdx( vec, val )
    for k = 1:length(vec)
        if( abs( vec(k) - val ) < 10^(-6) )
            idx=k;return;
        end
    end
    error('could not find the given val incide the vector vec!');
end