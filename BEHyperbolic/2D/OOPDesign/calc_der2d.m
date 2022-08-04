function [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d(vz,dvz,W,IDH,s_isdh,vdah,h,sx,alpha,beta,order)
 %================== 2nd
    VV = vdah;
    
    b = (alpha*beta*vz.*vz + (beta-1)*vz);
    
    %deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(vz,0,vdah,vdah);
    if(s_isdh == 3)  
        deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(vz,0,vdah,vdah);
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else
        deltab = W'* O4dh_gmI(b,0,vdah);    deltav  =  O4dh_gmI(vz,0,vdah);
        for j=1:sx
            VV(j,:) = pentsolv(IDH(j) - 1/12, [(1/12) (-16/12) IDH(j) (-16/12) 1/12],deltab(j,:));
        end
    end
    d2vz = (W*VV + deltav/h^2);
    if(order == 2) d3vz = vdah; d4vz = vdah; d5vz = vdah; d6vz = vdah; return; end
    
    %==================  3rd   
    
    b =(2*alpha*beta*dvz.*vz + (beta-1)*dvz);
    
    if(s_isdh == 3)  
        deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(dvz,0,vdah,vdah);
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else 
        deltab = W'* O4dh_gmI(b,0,vdah);    deltav  =  O4dh_gmI(dvz,0,vdah);
        for j=1:sx
            VV(j,:) = pentsolv(IDH(j) - 1/12, [(1/12) (-16/12) IDH(j) (-16/12) 1/12],deltab(j,:));
        end
    end
    d3vz = (W*VV + deltav/h^2);
    if(order == 3) d4vz = vdah; d5vz = vdah;  d6vz = vdah; return; end
    
    %==================  4th  
    
    b =( 2*alpha*beta*(dvz.*dvz + vz.*d2vz) + (beta-1)*d2vz);
    
    if(s_isdh == 3)  
        deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(d2vz,0,vdah,vdah);
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else 
        deltab = W'* O4dh_gmI(b,0,vdah);    deltav  =  O4dh_gmI(d2vz,0,vdah);
        for j=1:sx
            VV(j,:) = pentsolv(IDH(j) - 1/12, [(1/12) (-16/12) IDH(j) (-16/12) 1/12],deltab(j,:));
        end
    end
    d4vz = (W*VV + deltav/h^2);
    if(order == 4) d5vz = vdah;  d6vz = vdah; return; end
    
    %==================  5th
    
    b =( 2*alpha*beta*(3*dvz.*d2vz + vz.*d3vz) + (beta-1)*d3vz);

    if(s_isdh == 3)  
        deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(d3vz,0,vdah,vdah);
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else  
        deltab = W'* O4dh_gmI(b,0,vdah);    deltav  =  O4dh_gmI(d3vz,0,vdah);
        for j=1:sx
            VV(j,:) = pentsolv(IDH(j) - 1/12, [(1/12) (-16/12) IDH(j) (-16/12) 1/12],deltab(j,:));
        end
    end
    d5vz = (W*VV + deltav/h^2);
    if(order == 5)  d6vz = vdah; return; end
        %==================  6th
    
    b =( 2*alpha*beta*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) + (beta-1)*d4vz);

    if(s_isdh == 3)  
        deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(d4vz,0,vdah,vdah);
        for j=1:sx
            VV(j,:) = diag3solv([(-1) IDH(j) (-1)],deltab(j,:));
        end
    else  
        deltab = W'* O4dh_gmI(b,0,vdah);    deltav  =  O4dh_gmI(d4vz,0,vdah);
        for j=1:sx
            VV(j,:) = pentsolv(IDH(j) - 1/12, [(1/12) (-16/12) IDH(j) (-16/12) 1/12],deltab(j,:));
        end
    end
    d6vz = (W*VV + deltav/h^2);
    
end
    
function resm=dh_gmI(M,gm,vdah,vdah2)   %dh operator if gm = 0
      vdah(1:end-1,:) = M(2:end,:); 
      vdah(2:end,:) = vdah(2:end,:) + M(1:end-1,:);
      vdah2(:,1:end-1) = M(:,2:end);    
      vdah2(:,2:end) = vdah2(:,2:end) + M(:,1:end-1);
      
      resm = (vdah+vdah2) - (4+gm)*M;
end
function vdah=O4dh_gmI(M,gm,vdah)   %dh operator if gm = 0
    
      j=1;
      vdah(j,:) = - (2.5 -1/12 + gm)*M(j,:) + (16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      j=2;
      vdah(j,:) =    + (16/12)*M(j-1,:) - (2.5 + gm)*M(j,:) +  ( 16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      for j=3:size(M,1)-2
          vdah(j,:) = (-1/12)*M(j-2,:) + (16/12)*M(j-1,:) - (2.5 + gm)*M(j,:) +  (16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      end
      vdah(end-1,:) = (-1/12)*M(end-3,:)+ (16/12)*M(end-2,:) -  (2.5 + gm)*M(end-1,:) + (16/12)*M(end,:);
      vdah(end,:) =                        (-1/12)*M(end-2,:) + (16/12)*M(end-1,:) - (2.5 - 1/12 + gm)*M(end,:);
  
      
      j=1;
      vdah(:,j) = vdah(:,j) - (2.5 -1/12 )*M(:,j) + (16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      j=2;
      vdah(:,j) = vdah(:,j) +   (16/12)*M(:,j-1) - (2.5 )*M(:,j)  +  ( 16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      for j=3:size(M,2)-2
        vdah(:,j) = vdah(:,j) + (-1/12)*M(:,j-2) +  (16/12)*M(:,j-1) -   (2.5 )*M(:,j)  +  (16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      end
      vdah(:,end-1) =  vdah(:,end-1) + (-1/12)*M(:,end-3)+ (16/12)*M(:,end-2) - (2.5 )*M(:,end-1)  + (16/12)*M(:,end);
      vdah(:,end) =     vdah(:,end) +  (-1/12)*M(:,end-2) + (16/12)*M(:,end-1) - (2.5  - 1/12 )*M(:,end);
end  
      