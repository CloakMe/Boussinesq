clear;
    h=.1;
    x1=-20:h:20;
    v1=sech(x1);
    d4v1 = BMM ([-1 4 -6 4 -1],v1,-5);
    
    h=.05;
    x2=-20:h:20;
    v2=sech(x2);
    d4v2 = BMM ([-1 4 -6 4 -1],v2,-5);
    
    h=.025;
    x3=-20:h:20;
    v3=sech(x3);
    d4v3 = BMM ([-1 4 -6 4 -1],v3,-5);
    
    res1 = d4v1 - d4v2(1:2:end);
    res2 = d4v2(1:2:end) - d4v3(1:4:end);
    E1_1 = norm(res1(:),2)
    E2_1 = norm(res2(:),2)
plot(x1,d4v1,'k',x2,d4v2)
    conv = (log(E1_1)-log(E2_1))/log(2)