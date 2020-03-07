function T = exact_sol_p2(x,t,c,beta1,beta2,alfa)
    % exact solution for Double dispersion equation   with quadratic nonlinearity 
    k = -0.5 * sqrt((c*c-1)/(beta1*c*c-beta2));
    XK = k*(x-c*t);
    SK = sech(XK);
    SSK=SK.*SK;
    T=SSK*(1.5)*(c*c-1)/alfa;
    
end
