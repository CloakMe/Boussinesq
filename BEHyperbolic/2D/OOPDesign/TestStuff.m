figure(2)
mesh( this.x, this.y(1:5), vu(:,1:5)' );
title('vu')
xlabel('x');            ylabel('y');

pnts = 8;
domainUtilsP2 = BEDomainUtilsP2( this.x, this.y, this.order, this.beta, this.c, this.mu );
    newY = this.y(1) - pnts*this.h : this.h : this.y(1) - this.h ;
    sxx = length( this.x );
    syy = length( newY );
    ox1 = ones(1,syy);
    X = this.x'*ox1;
    if(sxx~=syy )
        oy1 = ones(1,sxx);
        Y = (newY'*oy1)';
    else
        Y = (newY'*ox1)';
    end
newLeft = domainUtilsP2.DersBnd( X, Y, this.tau + t(k) );

figure(4)

mesh( this.x, [ newY this.y(1:pnts)], [ newLeft(:,:,1) vu(:,1:pnts)]' );
title('vu-')
xlabel('x');            ylabel('y');