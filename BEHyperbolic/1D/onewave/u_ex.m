function [u]=u_ex(x,t,c,alpha,beta1,beta2)
%The formula below describes the exact solution for the following Boussinesq's Paradigm:
%d^2/dt^2(u) = /\(u) + beta1*(/\*d^2/dt^2)(u) - beta2*(/\^2)(u)+alpha*u^2.
%Here /\ stands for laplacian operator.


if (nargin == 5) beta2=0.5; end
if (nargin == 4) beta1=1.5; beta2=0.5;  end
if (nargin == 3) alpha=3; beta1=1.5; beta2=0.5; end

    beta=beta1/beta2;
    sx1 = size(x,1);
    sx2 = size(x,2);
    sx  = max(sx1,sx2);
    
if(sx > 1 && (sx1*sx2 == sx1 || sx1*sx2 == sx2))
    u=3*(c^2 - 1)*sech( sqrt( beta1*(c^2-1)/(beta1*c^2 - beta2) )*(x-c*t*sqrt(beta))/2 ).^2/(2*alpha);
else
    u=3*(c^2 - 1)*sech( sqrt( beta1*(c^2-1)/(beta1*c^2 - beta2) )*(x-c*t*sqrt(beta))/2 )^2/(2*alpha);
end