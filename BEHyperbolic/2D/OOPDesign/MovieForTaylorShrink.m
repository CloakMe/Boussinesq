mag = 120;
mag2=mag;

figure(20)
clear('F')
camva('manual');
for k = 1:sx
    if( -10^(-11)<x(k) )
        ijj=k;break;
    end
end
for k = 1:sy
   if( -10^(-11)<y(k))
       loo=k;break;
   end
end
tt = 1:0.1:20;    
for j = 1:size(tt,2)   
	
    xxx=x(ijj-mag+1:ijj+mag-1);
    yyy=y(loo-mag+1:loo+mag-1); 
    
    load(['SOL\MM_' num2str(tt(j)) '.mat']);
    [ crntSx, crntSy ] = size( vu );
    ij = ( crntSx + 1 )/2;
    lo = ( crntSy + 1 )/2;
    %vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
    vvu = vu(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag-1);
    mesh(xxx,yyy,vvu');
	%mesh(x,y,vu');
    title(['Solution at time: ',num2str(tt(j))]);
    xlabel('x'); ylabel('y');
    view(87,50);
    %view(1,90);
    axis([xxx(1) xxx(end) yyy(1) yyy(end) -.7, 2.5]);
    colorbar;
    caxis([-.7, 1.7]);
    F(:,j) = getframe;
    j=j+1;
    clear('v2u');clear('vu');
 end