 function Y = vesi_deltah(p,h)
Y=0*p;
Y(2:end-1)=(p(3:end)-2*p(2:end-1)+p(1:end-2))/(h*h);
end