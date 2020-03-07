function x=Comp(v)
    N=length(v);
    x=zeros(N/2,1);
    x(1)=v(1);
    x(end)=v(end);
    
    for i=2:(length(x)-1)
        x(i)=v(i*2-1);
    end
end