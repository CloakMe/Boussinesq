function MovieForBEHyperbolic( viewTypeX, viewTypeY, tt, x, y )
    mag = 120;
    %mag2=mag;
    h = x(2) - x(1);
    N = 200;
    M = 400;
    if(h > 0.15)
        N = 150;
        M = 300;
    end
    if(h > 0.35)
        N = 100;
        M = 200;
    end
    if(length(x) < 210 || length(y) < 420)
        x_st_idx = 1;
        x_st = x(1);
        x_end_idx = length(x);
        x_end = x(x_end_idx);

        y_st_idx = 1;
        y_st = y(y_st_idx);
        y_end_idx = length(y);
        y_end = y(y_end_idx);
    else
        x_st_idx = floor( (length(x) - N)/2 ) + 1;
        x_st = x(x_st_idx);
        x_end_idx = x_st_idx + N;
        x_end = x(x_end_idx);

        y_st_idx = floor( (length(y) - M)/2 ) + 1;
        y_st = y(y_st_idx);
        y_end_idx = y_st_idx + M;
        y_end = y(y_end_idx);
    end
        
    figure(6)
    clear('F')
    camva('manual');

    %tt = 1:0.1:20;    
    for j = 1:3:size(tt,2)   

        xxx = x(x_st_idx:x_end_idx);
        lenX = length( xxx );
        yyy = y(y_st_idx:y_end_idx); 
        lenY = length( yyy );
        load(['SOL\vu_' num2str(tt(j)) '.mat']);
        [ crntSx, crntSy ] = size( vu );
        %vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
        vvu = vu( x_st_idx:x_end_idx, y_st_idx :y_end_idx );
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