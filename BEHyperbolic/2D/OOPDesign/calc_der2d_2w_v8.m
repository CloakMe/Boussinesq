function [d2vz, d3vz, d4vz, d5vz, d6vz] =...
    calc_der2d_2w_v8(vz,dvz,W,IDH,socfd,sdah1,sdah2,vdah,h,t,sx,al,bt,c,c1,ord,sh,BNDR,M)

%================== 2nd
   VV = vdah;
   sy = size(vdah,2);
   [Ul] = drs_ex_bnd(BNDR.X.l,BNDR.Y.l,sdah1,t-sh,bt,c,c1,ord);%+drs_ex_bnd(Xl,Yl,sdah,t+sh,bt,c,c1,ord);
   [Ur] = drs_ex_bnd(BNDR.X.r,BNDR.Y.r,sdah1,t-sh,bt,c,c1,ord);
   [Ut] = drs_ex_bnd(BNDR.X.t,BNDR.Y.t,sdah2,t-sh,bt,c,c1,ord);
   [Ub] = drs_ex_bnd(BNDR.X.b,BNDR.Y.b,sdah2,t-sh,bt,c,c1,ord);
      
   [Nl] = NT_b(BNDR.X.l,BNDR.Y.l,sdah1,t-sh,al,bt,c,c1,ord);
   [Nr] = NT_b(BNDR.X.r,BNDR.Y.r,sdah1,t-sh,al,bt,c,c1,ord);
   [Nt] = NT_b(BNDR.X.t,BNDR.Y.t,sdah2,t-sh,al,bt,c,c1,ord);
   [Nb] = NT_b(BNDR.X.b,BNDR.Y.b,sdah2,t-sh,al,bt,c,c1,ord);
   
   [Bl] = u_ex_b(BNDR.X.l,BNDR.Y.l,sdah1,t-sh,bt,c,c1,ord);
   [Br] = u_ex_b(BNDR.X.r,BNDR.Y.r,sdah1,t-sh,bt,c,c1,ord);
   [Bt] = u_ex_b(BNDR.X.t,BNDR.Y.t,sdah2,t-sh,bt,c,c1,ord);
   [Bb] = u_ex_b(BNDR.X.b,BNDR.Y.b,sdah2,t-sh,bt,c,c1,ord);
   if(sh~=0)
       Ul = (Ul + drs_ex_bnd(BNDR.X.l,BNDR.Y.l,sdah1,t+sh,bt,c,c1,ord))*2/3;
       Ur = (Ur + drs_ex_bnd(BNDR.X.r,BNDR.Y.r,sdah1,t+sh,bt,c,c1,ord))*2/3;
       Ut = (Ut + drs_ex_bnd(BNDR.X.t,BNDR.Y.t,sdah2,t+sh,bt,c,c1,ord))*2/3;
       Ub = (Ub + drs_ex_bnd(BNDR.X.b,BNDR.Y.b,sdah2,t+sh,bt,c,c1,ord))*2/3;
       
       Nl = (Nl + NT_b(BNDR.X.l,BNDR.Y.l,sdah1,t+sh,al,bt,c,c1,ord))*2/3;
       Nr = (Nr + NT_b(BNDR.X.r,BNDR.Y.r,sdah1,t+sh,al,bt,c,c1,ord))*2/3;
       Nt = (Nt + NT_b(BNDR.X.t,BNDR.Y.t,sdah2,t+sh,al,bt,c,c1,ord))*2/3;
       Nb = (Nb + NT_b(BNDR.X.b,BNDR.Y.b,sdah2,t+sh,al,bt,c,c1,ord))*2/3;
   
       Bl = (Bl + u_ex_b(BNDR.X.l,BNDR.Y.l,sdah1,t+sh,bt,c,c1,ord))*2/3;
       Br = (Br + u_ex_b(BNDR.X.r,BNDR.Y.r,sdah1,t+sh,bt,c,c1,ord))*2/3;
       Bt = (Bt + u_ex_b(BNDR.X.t,BNDR.Y.t,sdah2,t+sh,bt,c,c1,ord))*2/3;
       Bb = (Bb + u_ex_b(BNDR.X.b,BNDR.Y.b,sdah2,t+sh,bt,c,c1,ord))*2/3;       
   end
   e = ones(sy,1);
   b = (al*bt*vz.*vz + (bt-1)*vz);
    
   %deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(vz,0,vdah,vdah);
   diagonal = 1:sy+1:sy^2;
   deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl(:,:,1),Nr(:,:,1),Nt(:,:,1),Nb(:,:,1),socfd) -...
   (BNDR.bndr_mat.l*Ut(:,:,1)' + BNDR.bndr_mat.r*Ub(:,:,1)' + (Ul(:,:,1)')*BNDR.bndr_mat.t +...
   (Ur(:,:,1)')*BNDR.bndr_mat.b)'); 

      deltav  =  Oh8_dh_gmI(vz,0,h,vdah,Bl(:,:,1),Br(:,:,1),Bt(:,:,1),Bb(:,:,1),socfd);
    for j=1:sx
        %M(1:sx+1:sx^2) = IDH(j)*e;
        M(diagonal) = IDH(j)*e;
        VV(j,:) = (M\deltab(j,:)');
    end
        
    d2vz = (W*VV + deltav/h^2);
    if(ord == 2) d3vz = vdah; d4vz = vdah; d5vz = vdah; d6vz = vdah; return; end
    
    %==================  3rd
    
    b =(2*al*bt*dvz.*vz + (bt-1)*dvz);
    
    deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl(:,:,2),Nr(:,:,2),Nt(:,:,2),Nb(:,:,2),socfd) -...
    (BNDR.bndr_mat.l*Ut(:,:,2)' + BNDR.bndr_mat.r*Ub(:,:,2)' + (Ul(:,:,2)')*BNDR.bndr_mat.t +...
    (Ur(:,:,2)')*BNDR.bndr_mat.b)'); 
    deltav  =  Oh8_dh_gmI(dvz,0,h,vdah,Bl(:,:,2),Br(:,:,2),Bt(:,:,2),Bb(:,:,2),socfd);

    for j=1:sx
        M(diagonal) = IDH(j)*e;
        VV(j,:) = (M\deltab(j,:)');
    end

    d3vz = (W*VV + deltav/h^2);
    if(ord == 3) d4vz = vdah; d5vz = vdah;  d6vz = vdah; return; end
    
    %==================  4th  
    
    b =( 2*al*bt*(dvz.*dvz + vz.*d2vz) + (bt-1)*d2vz);
    
    deltab = W'*( Oh8_dh_gmI(b,0,h,vdah,Nl(:,:,3),Nr(:,:,3),Nt(:,:,3),Nb(:,:,3),socfd) -...
    (BNDR.bndr_mat.l*Ut(:,:,3)' + BNDR.bndr_mat.r*Ub(:,:,3)' + (Ul(:,:,3)')*BNDR.bndr_mat.t +...
    (Ur(:,:,3)')*BNDR.bndr_mat.b)'); 
    deltav  =  Oh8_dh_gmI(d2vz,0,h,vdah,Bl(:,:,3),Br(:,:,3),Bt(:,:,3),Bb(:,:,3),socfd);

    for j=1:sx
        M(diagonal) = IDH(j)*e;
        VV(j,:) = (M\deltab(j,:)');
    end

    d4vz = (W*VV + deltav/h^2);
    if(ord == 4) d5vz = vdah;  d6vz = vdah; return; end
    
    %==================  5th
    
    b =( 2*al*bt*(3*dvz.*d2vz + vz.*d3vz) + (bt-1)*d3vz);

    deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl(:,:,4),Nr(:,:,4),Nt(:,:,4),Nb(:,:,4),socfd) -...
    (BNDR.bndr_mat.l*Ut(:,:,4)' + BNDR.bndr_mat.r*Ub(:,:,4)' + (Ul(:,:,4)')*BNDR.bndr_mat.t +...
    (Ur(:,:,4)')*BNDR.bndr_mat.b)'); 
    deltav  =  Oh8_dh_gmI(d3vz,0,h,vdah,Bl(:,:,4),Br(:,:,4),Bt(:,:,4),Bb(:,:,4),socfd);
        
    for j=1:sx
        M(diagonal) = IDH(j)*e;
        VV(j,:) = (M\deltab(j,:)');
    end

    d5vz = (W*VV + deltav/h^2);
    if(ord == 5)  d6vz = vdah; return; end
        %==================  6th
    
    b =( 2*al*bt*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) + (bt-1)*d4vz);

    deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl(:,:,5),Nr(:,:,5),Nt(:,:,5),Nb(:,:,5),socfd) -...
    (BNDR.bndr_mat.l*Ut(:,:,5)' + BNDR.bndr_mat.r*Ub(:,:,5)' + (Ul(:,:,5)')*BNDR.bndr_mat.t +...
    (Ur(:,:,5)')*BNDR.bndr_mat.b)'); 
    deltav  =  Oh8_dh_gmI(d4vz,0,h,vdah,Bl(:,:,5),Br(:,:,5),Bt(:,:,5),Bb(:,:,5),socfd);

    for j=1:sx
        M(diagonal) = IDH(j)*e;
        VV(j,:) = (M\deltab(j,:)');
    end

    d6vz = (W*VV + deltav/h^2);
    
end
    

function vdah=Oh8_dh_gmI(M,gm,h,vdah,uxst,uxend,uyst,uyend,socfd)   %dh operator if gm = 0
    vdah = Oh8_diff2y(M',gm,h,vdah',uxst',uxend',socfd)' + Oh8_diff2y(M,gm,h,vdah,uyst,uyend,socfd);
end

function vdah=Oh8_diff2x(M,gm,h,vdah,uxst,uxend,socfd)   %dh operator if gm = 0
    vdah= Oh8_diff2y(M',gm,h,vdah',uxst',uxend',socfd)';
end

function vdah=Oh8_diff2y(M,gm,h,vdah,uyst,uyend,socfd)   %dh operator if gm = 0

j=1;
vdah(:,j) = uyst(:,1:4)*socfd(1:4)' + M(:,j:j+4)*socfd(5:9)';
j=2;
vdah(:,j) = uyst(:,2:4)*socfd(1:3)' + M(:,j-1:j+4)*socfd(4:9)';
j=3;
vdah(:,j) = uyst(:,3:4)*socfd(1:2)' + M(:,j-2:j+4)*socfd(3:9)';
j=4;
vdah(:,j) =  uyst(:,4)*socfd(1) + M(:,j-3:j+4)*socfd(2:9)';
for j=5:size(M,2)-4
    vdah(:,j) = M(:,j-4:j+4)*socfd';
end

vdah(:,end-3) = M(:,end-7:end)*socfd(1:end-1)' + uyend(:,1)*socfd(end); 

vdah(:,end-2) = M(:,end-6:end)*socfd(1:end-2)' + uyend(:,1:2)*socfd(end-1:end)';
      
vdah(:,end-1) = M(:,end-5:end)* socfd(1:end-3)'+ uyend(:,1:3)*socfd(end-2:end)';

vdah(:,end) =  M(:,end-4:end)*socfd(1:end-4)' + uyend(:,1:4)*socfd(end-3:end)';

%{


cdf6=[1/90 	-3/20 	3/2 	-49/18 	3/2 	-3/20 	1/90];

cdf4=[-1/12 	4/3 	-5/2 	4/3 	-1/12 ];

cdf2=[-1 	2 	-1  ];

j=1;
vdah(:,j) = M(:,j:j+1)*cdf2(2:3)';
j=2;
vdah(:,j) =  M(:,j-1:j+2)*cdf4(2:5)';
j=3;
vdah(:,j) = M(:,j-2:j+3)*cdf6(2:7)';
j=4;
vdah(:,j) =M(:,j-3:j+4)*socfd(2:9)';
for j=5:size(M,2)-4
    vdah(:,j) = M(:,j-4:j+4)*socfd';
end

vdah(:,end-3) = M(:,end-7:end)*socfd(1:end-1)' ; 

vdah(:,end-2) = M(:,end-5:end)*cdf6(1:end-1)' ;
      
vdah(:,end-1) = M(:,end-3:end)* cdf4(1:end-1)';

vdah(:,end) =  M(:,end-1:end)*cdf2(1:2)';
%}
end