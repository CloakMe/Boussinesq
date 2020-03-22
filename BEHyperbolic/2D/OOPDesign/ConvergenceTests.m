clear;
bndCutSizeX = 0;
bndCutSizeY = 0;

yo(1,:) = '20';
yo(2,:) = '10';
yo(3,:) = '05';
%_tau05
yo = cellstr(yo);
    
for jl = 1:3

    cellStr = strcat('SavedWorkspaces\ZeroBoundary\EnergySave_40_ZB1_bt3_c052_h0', yo(jl), '_O(h^2)' );
    %cellStr = strcat('Hyp_40_bt1_c090_h0', yo(jl), '_O(h^4)' );
    warning('off','all');
    load (  cellStr{1} );
    warning('on','all');
    %fprintf('tau = %f\n', tau);
    %fprintf( 'h = 0.%s \n',  yo{jl});
    bndPtsRemX = 0.0;
    bndPtsRemY = 0.0;
    if(jl==3)
        bndPtsRemX = bndCutSizeX/( hc/2.0 );
        bndPtsRemY = bndCutSizeY/( hc/2.0 );
    else
        bndPtsRemX = bndCutSizeX/waveFactory.h;
        bndPtsRemY = bndCutSizeY/waveFactory.h;   
    end
    %sum(tauVector)
    %xx = waveFactory.x( bndPtsRemX+1:end-bndPtsRemX );
    %yy = waveFactory.y( bndPtsRemY+1:end-bndPtsRemY );
    
    vl = vl( bndPtsRemX+1:end-bndPtsRemX, bndPtsRemY+1:end-bndPtsRemY );
    if(jl==1)
       vb = vl(bndPtsRemX+1:end-bndPtsRemX,bndPtsRemX+1:end-bndPtsRemX); 
       hb=waveFactory.h
       taub = tau
       ENb = EN;
    end
    if(jl==2)
       vc = vl(bndPtsRemX+1:end-bndPtsRemX,bndPtsRemX+1:end-bndPtsRemX); 
       hc=waveFactory.h
       tauc = tau
       ENc = EN;
    end
    if(jl==3)
       vd = vl(bndPtsRemX+1:end-bndPtsRemX,bndPtsRemX+1:end-bndPtsRemX);  
       taud = tau
       ENd = EN;
    end
    clear('U');clear('cellStr');clear('waveFactory');
end

clear('jl');
res_1 = vb - vc(1:2:end,1:2:end);
res_2 = vc - vd(1:2:end,1:2:end);

sizeSmall = size( ENb, 2 )
sizeMedium = size( ENc, 2 )
sizeBig = size( ENd, 2 )
if( sizeSmall ~= 1 && 2*sizeSmall+1 == sizeMedium && ...
    sizeMedium ~= 1 && 2 * sizeMedium + 1 == sizeBig )
    ENb = [ 0, ENb ];
    ENc = [ 0, ENc ];
    ENd = [ 0, ENd ];
    fprintf('adjusting energy sizes\n');
end

EN_1 = ENb - ENc(1:2:end);
EN_2 = ENc - ENd(1:2:end);


norm0402_L2 = hb*norm(res_1(:),2);
norm0201_L2 = hc*norm(res_2(:),2);
fprintf('Solution Convergence:\n');

fprintf('||v_04 - v_02||_L2 = %.6f \n', norm0402_L2);
fprintf('||v_02 - v_01||_L2 = %.6f \n', norm0201_L2); %%.8e
conv_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2);
fprintf('Conv_L2 = %.8e \n\n', conv_L2);

norm0402_Inf = max(max(abs(res_1(:))));
norm0201_Inf = max(max(abs(res_2(:))));
fprintf('||v_04 - v_02||_Inf = %.6f \n', norm0402_Inf);
fprintf('||v_02 - v_01||_Inf = %.6f \n', norm0201_Inf);
conv_Inf = log(abs(norm0402_Inf/norm0201_Inf))/log(2);
fprintf('Conv_Inf = %.8e \n', conv_Inf);

fprintf('\nEnergy Convergence:\n');

norm0402_L2 = taub*norm(EN_1(:),2);
norm0201_L2 = tauc*norm(EN_2(:),2);
fprintf('||EN_04 - EN_02||_L2 = %.6f \n', norm0402_L2);
fprintf('||EN_02 - EN_01||_L2 = %.6f \n', norm0201_L2); %%.8e
conv_L2 = log(abs(norm0402_L2/norm0201_L2))/log(2);
fprintf('Conv_L2 = %.8e \n\n', conv_L2);

norm0402_Inf = max(max(abs(EN_1)));
norm0201_Inf = max(max(abs(EN_2)));
fprintf('||EN_04 - EN_02||_Inf = %.6f \n', norm0402_Inf);
fprintf('||EN_02 - EN_01||_Inf = %.6f \n', norm0201_Inf); %%.8e
conv_L2 = log(abs(norm0402_Inf/norm0201_Inf))/log(2);
fprintf('Conv_Inf = %.8e \n\n', conv_L2);

clear;