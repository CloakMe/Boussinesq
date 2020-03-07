function [X] = deltahat(p,h)
N=length(p);
X=0*p;  
znam=4*h*h;
X(1)=(p(3)-2*p(1)+p(N-1))/znam;
X(2)=(p(4)-2*p(2)+p(N))/znam;
X(3:N-2)=(p(5:N)-2*p(3:N-2)+p(1:N-4))/znam;
X(N-1)=(p(1)-2*p(N-1)+p(N-3))/znam;
X(N)=(p(2)-2*p(N)+p(N-2))/znam;
end
