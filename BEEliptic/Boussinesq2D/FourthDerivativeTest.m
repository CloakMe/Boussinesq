function [  ] = FourthDerivativeTest(  )
    
    %load (['..\SavedSolutions\who22_bt3_c0' num2str(gg) '_h02']);  %gg = [1 3 5]
    %testConvIterMethod
    load (['TestConvIterMethod\who_30_bt3_c010_h02_ord2']);
    u_xxxx04 =XDer(XDer(bigU,derivative.second),derivative.second);
    u_xxyy04 =XDer(YDer(bigU,derivative.second),derivative.second);
    u_yyyy04 =YDer(YDer(bigU,derivative.second),derivative.second);
    h1=h;
    
    clearvars -except u_xxxx04 u_xxyy04 u_yyyy04 h1;
    
    load (['TestConvIterMethod\who_30_bt3_c010_h01_ord2']);
    u_xxxx02 =XDer(XDer(bigU,derivative.second),derivative.second);
    u_xxyy02 =XDer(YDer(bigU,derivative.second),derivative.second);
    u_yyyy02 =YDer(YDer(bigU,derivative.second),derivative.second);
    h2=h;
        
    clearvars -except u_xxxx04 u_xxyy04 u_yyyy04 u_xxxx02 u_xxyy02 u_yyyy02 h1 h2;
    
    load (['TestConvIterMethod\who_30_bt3_c010_h005_ord2']);
    u_xxxx01 =XDer(XDer(bigU,derivative.second),derivative.second);
    u_xxyy01 =XDer(YDer(bigU,derivative.second),derivative.second);
    u_yyyy01 =YDer(YDer(bigU,derivative.second),derivative.second);
    
    %=========================================================
    res_xxxx0402 = u_xxxx04(13:end-12,13:end-12) - u_xxxx02(25:2:end-24,25:2:end-24);
    res_xxxx0201 = u_xxxx02(13:end-12,13:end-12) - u_xxxx01(25:2:end-24,25:2:end-24);

    norm0402_L2 = h1*norm(res_xxxx0402(:),2);
    norm0201_L2 = h2*norm(res_xxxx0201(:),2);
    fprintf('||u_xxxx_04 - u_xxxx_02||_L2 = %.4e \n', norm0402_L2);
    fprintf('||u_xxxx_02 - u_xxxx_01||_L2 = %.4e \n', norm0201_L2);
    conv_xxxx_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2)

    norm0402_Inf = max(max(abs(res_xxxx0402(:))));
    norm0201_Inf = max(max(abs(res_xxxx0201(:))));
    fprintf('||u_xxxx_04 - u_xxxx_02||_Inf = %.4e \n', norm0402_Inf);
    fprintf('||u_xxxx_02 - u_xxxx_01||_Inf = %.4e \n', norm0201_Inf);
    conv_xxxx_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2)
    
    %=========================================================
    res_xxyy0402 = u_xxyy04(13:end-12,13:end-12) - u_xxyy02(25:2:end-24,25:2:end-24);
    res_xxyy0201 = u_xxyy02(13:end-12,13:end-12) - u_xxyy01(25:2:end-24,25:2:end-24);

    norm0402_L2 = h1*norm(res_xxyy0402(:),2);
    norm0201_L2 = h2*norm(res_xxyy0201(:),2);
    fprintf('||u_xxyy_04 - u_xxyy_02||_L2 = %.4e \n', norm0402_L2);
    fprintf('||u_xxyy_02 - u_xxyy_01||_L2 = %.4e \n', norm0201_L2);
    conv_xxyy_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2)

    norm0402_Inf = max(max(abs(res_xxyy0402(:))));
    norm0201_Inf = max(max(abs(res_xxyy0201(:))));
    fprintf('||u_xxyy_04 - u_xxyy_02||_Inf = %.4e \n', norm0402_Inf);
    fprintf('||u_xxyy_02 - u_xxyy_01||_Inf = %.4e \n', norm0201_Inf);
    conv_xxyy_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2)
    
    %=========================================================
  
    res_yyyy0402 = u_yyyy04(13:end-12,13:end-12) - u_yyyy02(25:2:end-24,25:2:end-24);
    res_yyyy0201 = u_yyyy02(13:end-12,13:end-12) - u_yyyy01(25:2:end-24,25:2:end-24);

    norm0402_L2 = h1*norm(res_yyyy0402(:),2);
    norm0201_L2 = h2*norm(res_yyyy0201(:),2);
    fprintf('||u_yyyy04 - u_yyyy02||_L2 = %.4e \n', norm0402_L2);
    fprintf('||u_yyyy02 - u_yyyy01||_L2 = %.4e \n', norm0201_L2);
    conv_yyyy_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2)

    norm0402_Inf = max(max(abs(res_yyyy0402(:))));
    norm0201_Inf = max(max(abs(res_yyyy0201(:))));
    fprintf('||u_yyyy04 - u_yyyy02||_Inf = %.4e \n', norm0402_Inf);
    fprintf('||u_yyyy02 - u_yyyy01||_Inf = %.4e \n', norm0201_Inf);
    conv_yyyy_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2)
    
    
        %fprintf('u_xxxx(0,0) = %.4e \n', u_xxxx(zeroX,zeroY));
        %fprintf('u_xxyy(0,0) = %.4e \n', u_xxyy(zeroX,zeroY));
        %fprintf('u_yyyy(0,0) = %.4e \n\n', u_yyyy(zeroX,zeroY));

    
end

