clear;
%Eliptic
bndCutSizeX = 0;
bndCutSizeY = 0;
bndDirectory = 'ZeroBoundary\'; % 'ZeroBoundary\' 'WithBoundary\'
ICType = 'ChristovIC_'; % 'ChristovIC_' 'Natali_'
domLength = '40_'; %'80_' '40_' '82_'
paramString = 'bt3_c052'; %'bt3_c052' 'bt1_c090'
orderString = '2';
%============================================
fprintf(paramString);
fprintf('!\n');
yo(1,:) = '20';
yo(2,:) = '10';
yo(3,:) = '05';
if( strcmp(paramString, 'bt1_c090') == 1 )
    yo(1,:) = '40';
    yo(2,:) = '20';
    yo(3,:) = '10';
end
yo = cellstr(yo);
bndString = 'ZB0_';
if( strcmp(bndDirectory, 'ZeroBoundary\') == 1 )
    bndString = 'ZB1_';
end
if( strcmp(ICType, 'EnergySave_') == 1 && strcmp(orderString, '2') == 0 )
    fprintf('EnergySave solution; setting order = 2!\n');
    orderString = '2';
end
directory = strcat(bndDirectory, ICType, domLength, paramString, '\', 'Oh', orderString, '\');

for jl = 1:3            
    strName = strcat(directory, ICType, domLength, bndString, paramString, '_h0', yo(jl), '_O(h^', orderString, ').mat' );
    load (  strName{1}, 'U', 'h' );
    %sum(tauVector)
    bndPtsRemX = 0;
    bndPtsRemY = 0;
    bndPtsRemX = bndCutSizeX/h;
    bndPtsRemY = bndCutSizeY/h;
    if(jl==1)
       vb = U(1:end-bndPtsRemX,1:end-bndPtsRemX); hb=h; 
    end
    if(jl==2)
       vc = U(1:end-bndPtsRemX,1:end-bndPtsRemX); hc=h;
    end
    if(jl==3)
       vd = U(1:end-bndPtsRemX,1:end-bndPtsRemX); 
    end
    %PlotAssymptotics(x,y,h,zeroX,zeroY,bigU,0);
    clear('U');clear('h');clear('strName');
end

clear('jl');
res_1 = vb - vc(1:2:end,1:2:end);
res_2 = vc - vd(1:2:end,1:2:end);

norm0402_L2 = hb*norm(res_1(:),2);
norm0201_L2 = hc*norm(res_2(:),2);
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

clear;