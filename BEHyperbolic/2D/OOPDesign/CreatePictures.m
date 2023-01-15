function CreatePictures( viewTypeX, viewTypeY, tt, x, y )
    mag = 120;
    %mag2=mag;
    h = x(2) - x(1);
    N = 375;
    M = 375;
    if(h > 0.15)
        N = 150;
        M = 400;
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
    %j = 1:3:size(tt,2)  
    for j = 0:4:20   

        xxx = x(x_st_idx:x_end_idx);
        lenX = length( xxx );
        yyy = y(y_st_idx:y_end_idx); 
        lenY = length( yyy );
        load(['SOL_bt3_c030\vz_' num2str(j) '.mat']);
        [ crntSx, crntSy ] = size( vz );
        %vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
        vvz = vz( x_st_idx:x_end_idx, y_st_idx :y_end_idx );
        surf(xxx,yyy,vvz');
        shading interp
        %mesh(x,y,vu');
        title(['t =  ',num2str(j)]);
        xlabel('t','FontSize',18); % ylabel('FontSize',18);
        set(gca,'FontSize',18);
        %view(1,90);
        %axis([xxx(1) xxx(end) yyy(1) yyy(end) -.7, 2.5]);
        %colorbar;
        axis tight;
        %caxis([-.7, 1.7]);
        view( viewTypeX, viewTypeY );
        
        %F(:,j) = getframe;
        %j=j+1;
        clear('v2u');clear('vz');        
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