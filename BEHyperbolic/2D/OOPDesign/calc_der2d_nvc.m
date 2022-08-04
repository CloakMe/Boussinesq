function [d2vz, d3vz, d4vz, d5vz] = calc_der2d_nvc(vz,dvz,W,IDH,s_isdh,vdah,h,sx,al,beta1,beta2,order)
    %================== 2nd time der
    VV = vdah;
    
    b=(vz + al*vz.*vz - (beta2/h^2)*dh_gmI(vz,0,vdah,vdah))/beta1;
        
    deltab = W'* dh_gmI(b,0,vdah,vdah);
            
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else 
        %FPP for five diag matrix
    end

    d2vz = W*VV;
    if(order == 2) d3vz = vdah; d4vz = vdah; d5vz = vdah; return; end
    
    %==================  3rd time der   
    
    b=(dvz + 2*al*dvz.*vz - (beta2/h^2)*dh_gmI(dvz,0,vdah,vdah))/beta1;
    
    deltab = W'* dh_gmI(b,0,vdah,vdah);    
    
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else 
        %FPP for five diag matrix
    end
    d3vz = W*VV;
    if(order == 3) d4vz = vdah; d5vz = vdah; return; end
    
    %================== 4th time der
    
    b=(d2vz + 2*al*(dvz.*dvz + vz.*d2vz) - (beta2/h^2)*dh_gmI(d2vz,0,vdah,vdah))/beta1;
    
    deltab = W'* dh_gmI(b,0,vdah,vdah);
    
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else  
    end
    d4vz = W*VV;
    if(order == 4) d5vz = vdah; return; end
    
    %==================  5th time der
    
    b=(d3vz + 2*al*(3*dvz.*d2vz + vz.*d3vz) - (beta2/h^2)*dh_gmI(d3vz,0,vdah,vdah))/beta1;

    deltab = W'* dh_gmI(b,0,vdah,vdah); 
    
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else  
    end
    d5vz = W*VV;
    
    
function resm=dh_gmI(M,gm,vdah,vdah2)   %dh operator if gm = 0
      vdah(1:end-1,:) = M(2:end,:); 
      vdah(2:end,:) = vdah(2:end,:) + M(1:end-1,:);
      vdah2(:,1:end-1) = M(:,2:end);    
      vdah2(:,2:end) = vdah2(:,2:end) + M(:,1:end-1);
      
      resm = (vdah+vdah2) - (4+gm)*M;