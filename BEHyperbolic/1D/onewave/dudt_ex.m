function [du]=dudt_ex(x,t,c,alpha,beta1,beta2)
%The formula below describes the time derivative of a function "u"
%where "u" is the exact solution for the following Boussinesq's Paradigm:
%d^2/dt^2(u) = /\(u) + beta1*(/\*d^2/dt^2)(u) - beta2*(/\^2)(u)+alpha*u^2.
%Here /\ stands for laplacian operator.

if (nargin == 5) beta2=0.5; end
if (nargin == 4) beta1=1.5; beta2=0.5;  end
if (nargin == 3) alpha=3; beta1=1.5; beta2=0.5; end

    beta=beta1/beta2;
    
   du= 3*(c^2-1)*( sech( sqrt( beta1*(c^2-1)/(beta1*c^2-beta2) ) * ( x-c*t*sqrt(beta) ) / 2 ) ).^2 .*...
  tanh( sqrt( beta1*(c^2-1)/(beta1*c^2-beta2) ) * ( x-c*t*sqrt(beta) ) /2 ) *sqrt( beta1*(c^2-1)/(beta1*c^2-beta2) ) * c * sqrt(beta) /(2*alpha);

   