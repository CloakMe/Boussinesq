function Y=B(u,alfa,beta2,h)
% samo pri kvadratichna nelinejnost
nonlinear=u.*u;
Y= xhat(u,h)-beta2*xhatxhatxhat(u,h)+alfa*xhat(nonlinear,h);
end