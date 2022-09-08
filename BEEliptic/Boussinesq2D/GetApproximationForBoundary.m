function [boundaryU2, approxBoundaryF_xx, approxBoundaryF_yy] = GetApproximationForBoundary(x,y,h,c)

   [X,Y]=AugDomain(x,y,h);
   c12 = 1-c^2;
   %boundaryU1 = 2 * (sqrt(c12) * X .* Y)./(c12*X.^2 + Y.^2).^2;    
   boundaryU2 = (c12 * X.^2 - Y.^2)./(c12*X.^2 + Y.^2).^2; 
   
   approxBoundaryF_xx =  6*c12* (Y.^4 - 6*c12*X.^2.*Y.^2 + c12 * X.^4) ./ (Y.^2 + c12*X.^2).^4;
   approxBoundaryF_yy = -6 *    (Y.^4 - 6*c12*X.^2.*Y.^2 + c12 * X.^4) ./ (Y.^2 + c12*X.^2).^4;
   %boundaryU=(c12^2* X.^4 - 6*c12 * X.^2 .* Y.^2 + Y.^4)./(c12*X.^2+Y.^2).^4;  
end
%extending the domain by adding four points top and right 
function [X,Y]=AugDomain(x,y,h)
    
    x2 = [x x(end)+h x(end)+2*h x(end)+3*h x(end)+4*h];
    y2 = [y y(end)+h y(end)+2*h y(end)+3*h y(end)+4*h];
   
    sxn = length(x2);
    syn = length(y2);
    ox1 = ones(1,syn);
    X = x2'*ox1;
   if(sxn~=syn)
        oy1 = ones(1,sxn);
        Y = (y2'*oy1)';
    else
        Y = (y2'*ox1)';
   end
end