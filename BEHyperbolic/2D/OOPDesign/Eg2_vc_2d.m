function [e] = Eg2_vc_2d(vmo,vz,vpo,IDHo,W,h,tau,al,bt,sx,sy,s_idh,vdah,VV,vec1)
   vt = (vpo - vz)/tau;
   wvt = W'*vt;
   idhv = -dh_gmI_v2(vz+vpo,h^2,s_idh,vdah,vdah)/h^2; %/h^2
   
  if(s_idh == 3)  
        for j=1:sx
            VV(j,:) = BEUtilities.TridiagSolv([(-1) IDHo(j) (-1)]/h^2,wvt(j,:));
        end
    else
        for j=1:sx
            VV(j,:) = BEUtilities.PentSolv(IDHo(j), [(1/12) (-16/12) IDHo(j) (-16/12) 1/12]/h^2,wvt(j,:));
        end
  end
    vec1 = W*VV;  %/h^2
    
    Le= h^2*(sum( sum(vec1.*vt) ))  + h^2*sum(sum(vt.*vt))  +  (h^2/4)*sum(sum(idhv.*(vz+vpo)));
       
    NLe = sum(sum((al*bt/3)*vz.^3 + ((bt-1)/2)*vz.^2))*(sx*sy*h^2)/(sx*sy);
    NLe =  (h^2*al*bt/3)*(  sum(sum( vz.^3 ))+ sum(sum( vpo.^3 ))  ) + (h^2*(bt-1)/2)*(  sum(sum(vz.^2)) + sum(sum(vpo.^2))  ); 
       
    e = Le + NLe;
  
end

function vdah=dh_gmI_v2(M,gm,s_idh,vdah,vdah2)   %dh operator if gm = 0
    if(s_idh == 3)
      vdah(1:end-1,:) = M(2:end,:); 
      vdah(2:end,:) = vdah(2:end,:) + M(1:end-1,:);
      vdah2(:,1:end-1) = M(:,2:end);    
      vdah2(:,2:end) = vdah2(:,2:end) + M(:,1:end-1);
      
      vdah = (vdah+vdah2) - (4+gm)*M;
    else
      
      j=1;
      vdah(j,:) = - (2.5 + gm)*M(j,:) + (16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      j=2;
      vdah(j,:) =    + (16/12)*M(j-1,:) - (2.5 + gm)*M(j,:) +  ( 16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      for j=3:size(M,1)-2
          vdah(j,:) = (-1/12)*M(j-2,:) + (16/12)*M(j-1,:) - (2.5 + gm)*M(j,:) +  (16/12)*M(j+1,:) + (-1/12)*M(j+2,:);
      end
      vdah(end-1,:) = (-1/12)*M(end-3,:)+ (16/12)*M(end-2,:) -  (2.5 + gm)*M(end-1,:) + (16/12)*M(end,:);
      vdah(end,:) =                        (-1/12)*M(end-2,:) + (16/12)*M(end-1,:) - (2.5 + gm)*M(end,:);
        
      j=1;
      vdah(:,j) = vdah(:,j) - (2.5 )*M(:,j) + (16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      j=2;
      vdah(:,j) = vdah(:,j) +   (16/12)*M(:,j-1) - (2.5 )*M(:,j)  +  ( 16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      for j=3:size(M,2)-2
        vdah(:,j) = vdah(:,j) + (-1/12)*M(:,j-2) +  (16/12)*M(:,j-1) -   (2.5 )*M(:,j)  +  (16/12)*M(:,j+1) + (-1/12)*M(:,j+2);
      end
      vdah(:,end-1) =  vdah(:,end-1) + (-1/12)*M(:,end-3)+ (16/12)*M(:,end-2) - (2.5 )*M(:,end-1)  + (16/12)*M(:,end);
      vdah(:,end) =     vdah(:,end) +  (-1/12)*M(:,end-2) + (16/12)*M(:,end-1) - (2.5 )*M(:,end);
      
    end 
end
