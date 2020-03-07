function X=xhat(p,h)
N=length(p);
X=0*p;
X(1)=(p(2)-p(N))/(2*h);
X(2:N-1)=(p(3:N)-p(1:N-2))/(2*h);
X(N)=(p(1)-p(N-1))/(2*h);
end