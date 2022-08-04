function [d2vz, d3vz, d4vz, d5vz] = calc_der2d_v2(vz,dvz,W,IDH,s_isdh,vdah,h,sx,al,beta1,beta2,order)
 %==================
    VV = vdah;
    
    b=(al*vz.*vz + (beta2/h^2)*(-dh_gmI(vz,h^2/beta2,vdah)));
    
    deltab = W'* dh_gmI(b,0,vdah)/beta1; 
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else 
    end
    d2vz = W*VV ;
    if(order == 2) d3vz = vdah; d4vz = vdah; d5vz = vdah; return; end
    
    %==================    
    
    b=( 2*al*dvz.*vz + (beta2/h^2)*(-dh_gmI(dvz,h^2/beta2,vdah)));
    deltab = W'* dh_gmI(b,0,vdah)/beta1;    
    
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else 
    end
    d3vz = W*VV ;
    if(order == 3) d4vz = vdah; d5vz = vdah; return; end
    
    %==================    
    
    b=(2*al*(dvz.*dvz + vz.*d2vz) + (beta2/h^2)*(-dh_gmI(d2vz,h^2/beta2,vdah)));
    deltab = W'* dh_gmI(b,0,vdah)/beta1;    
    
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else  
    end
    d4vz = W*VV ;
    if(order == 4) d5vz = vdah; return; end
    
    %==================  
    
    b=(2*al*(3*dvz.*d2vz + vz.*d3vz) + (beta2/h^2)*(-dh_gmI(d3vz, h^2/beta2 ,vdah)));
    deltab = W'* dh_gmI(b,0,vdah)/beta1;    
    
    if(s_isdh == 3)  
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else  
    end
    d5vz = W*VV ;

    %==================  
    
function vdah=dh_gmI(M,gm,vdah)   %dh operator if gm = 0

      j=1;
      vdah(j,:) = - (2+gm)*M(j,:) +  M(j+1,:);
      for j=2:size(M,1)-1
          vdah(j,:) =  M(j-1,:) - (2+gm)*M(j,:) +  M(j+1,:);
      end
      j=size(M,1);
      vdah(j,:) =  M(j-1,:) - (2+gm)*M(j,:) ;
      
      
      j=1;
      vdah(:,j) = vdah(:,j) - (2)*M(:,j) +  M(:,j+1);
      for j=2:size(M,1)-1
          vdah(:,j) = vdah(:,j)  +  M(:,j-1) - (2)*M(:,j) +M(:,j+1);
      end
      j=size(M,1);
      vdah(:,j) =  vdah(:,j) +M(:,j-1) - (2)*M(:,j) ;