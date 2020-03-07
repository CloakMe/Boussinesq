function [X]=xhatxhatxhat(p,h)

    N=length(p);
    X=0*p;
    znam=8*h*h*h;
    X(1)=(p(4)-3*p(2)+3*p(N)-p(N-2))/znam;
    X(2)=(p(5)-3*p(3)+3*p(1)-p(N-1))/znam;
    X(3)=(p(6)-3*p(4)+3*p(2)-p(N))/znam;
    X(4:N-3)=(p(7:N)-3*p(5:N-2)+3*p(3:N-4)-p(1:N-6))/znam;
    X(N-2)=(p(1)-3*p(N-1)+3*p(N-3)-p(N-5))/znam;
    X(N-1)=(p(2)-3*p(N)+3*p(N-2)-p(N-4))/znam;
    X(N)=(p(3)-3*p(1)+3*p(N-1)-p(N-3))/znam;
    
end