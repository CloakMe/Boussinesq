domutils = BEDomainUtilsP2( x, y, 4, 3, .45, waveFactory.mu );
mag2=mag;
figure(6)
clear('F')
camva('manual');
for j = 1:size(tt,2)
    if(ver == 6)
         load(['SOL\MM2_' num2str(tt(j)) '.mat']);
         %vv1u = v1u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
         vv2u = v2u(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1);
         mesh(xx,yy,vv2u');
         %mesh(x(1:xstep:end),y(1:xstep:end),v2u(1:xstep:end,1:xstep:end)');

    else             
        load(['SOL\MM_' num2str(tt(j)) '.mat']);
        vvu = vu(ij-mag+1:ij+mag-1,lo-mag+1:lo+mag2-1); 
        figure( 6 );
        mesh(x,y,vu');
        title(['Solution at time: ',num2str(tt(j))]);
        xlabel('x'); ylabel('y');
        axis([x(1) x(end) y(1) y(end) -.001, 0.001]);
        caxis([-.001, 0.001]);
        %view(87,50);
        view(90,90);
        %mesh(x(1:xstep:end),y(1:xstep:end),vu(1:xstep:end,1:xstep:end)');
        figure( 7 );
        dersBnd = domutils.GetDersBndMid( tt(j) );
        mesh( x, y, dersBnd(:,:,1)' );
        title(['Bnd Fun at time: ',num2str(tt(j))]);
        xlabel('x'); ylabel('y');
        axis([x(1) x(end) y(1) y(end) -.001, 0.001]);
        caxis([-.001, 0.001]);
        view(90,90);

    end
    %mesh(x,y,vu');
    F(:,j) = getframe;
    j=j+1;
    clear('v2u');clear('vu');
 end