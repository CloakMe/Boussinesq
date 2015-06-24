figure(2)
mesh(this.x,this.y,deltaTimeder')
title('vu , t=this.tau')
xlabel('x');            ylabel('y');
xInd = 1:10;
endy = length(this.y);
yInd = endy-9:endy;
figure(2)
mesh(this.x(xInd),this.y(yInd),deltaTimeder(xInd,yInd)')
title('vu , t=this.tau')
xlabel('x');            ylabel('y');

figure(3)
mesh(this.x(xInd),this.y(yInd),vu(xInd,yInd)')
title('vu , t=this.tau')
xlabel('x');            ylabel('y');


figure(2)
        mesh( X(xInd,1)', Y(1,yInd), (extrapolateEdge(xInd,yInd) - edge(xInd,yInd))' );
        figure(3)
        mesh( X(xInd,1)', Y(1,yInd), (extrapolateEdge(xInd,yInd))' );
        figure(4)
        mesh( X(xInd,1)', Y(1,yInd), ( edge(xInd,yInd))' );
        
        figure(5)
        mesh( this.x(xInd)', this.x(yInd), ( vz(xInd,yInd))' );
        figure(6)
        mesh( this.x(xInd)', this.x(yInd), ( vz(xInd,yInd))' );