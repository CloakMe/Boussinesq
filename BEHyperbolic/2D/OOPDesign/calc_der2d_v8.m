function [d2vz, d3vz, d4vz, d5vz, d6vz] = calc_der2d_v8(vz,dvz,W,IDH,socfd,vdah,h,t,sx,al,bt,c,c1,ord,Xl,Xr,Xt,Xb,Yl,Yr,Yt,Yb,DRs_l,DRs_r,DRs_t,DRs_b,M)

%================== 2nd
   VV = vdah;
   sy = size(vdah,2);
   [Ul2, Ul3, Ul4, Ul5, Ul6, Ul7, Ul8] = drs_ex_bnd(Xl,Yl,t,bt,c,c1,ord);
   [Ur2, Ur3, Ur4, Ur5, Ur6, Ur7, Ur8] = drs_ex_bnd(Xr,Yr,t,bt,c,c1,ord);
   [Ut2, Ut3, Ut4, Ut5, Ut6, Ut7, Ut8] = drs_ex_bnd(Xt,Yt,t,bt,c,c1,ord);
   [Ub2, Ub3, Ub4, Ub5, Ub6, Ub7, Ub8] = drs_ex_bnd(Xb,Yb,t,bt,c,c1,ord);
   
   [Nl2, Nl3, Nl4, Nl5, Nl6, Nl7, Nl8] = NT_b(Xl,Yl,t,al,bt,c,c1,ord);
   [Nr2, Nr3, Nr4, Nr5, Nr6, Nr7, Nr8] = NT_b(Xr,Yr,t,al,bt,c,c1,ord);
   [Nt2, Nt3, Nt4, Nt5, Nt6, Nt7, Nt8] = NT_b(Xt,Yt,t,al,bt,c,c1,ord);
   [Nb2, Nb3, Nb4, Nb5, Nb6, Nb7, Nb8] = NT_b(Xb,Yb,t,al,bt,c,c1,ord);
   
   [Bl2, Bl3, Bl4, Bl5, Bl6, Bl7, Bl8] = u_ex_b(Xl,Yl,t,bt,c,c1,ord);
   [Br2, Br3, Br4, Br5, Br6, Br7, Br8] = u_ex_b(Xr,Yr,t,bt,c,c1,ord);
   [Bt2, Bt3, Bt4, Bt5, Bt6, Bt7, Bt8] = u_ex_b(Xt,Yt,t,bt,c,c1,ord);
   [Bb2, Bb3, Bb4, Bb5, Bb6, Bb7, Bb8] = u_ex_b(Xb,Yb,t,bt,c,c1,ord);
   e = ones(sy,1);
   b = (al*bt*vz.*vz + (bt-1)*vz);
    
    %deltab = W'* dh_gmI(b,0,vdah,vdah);    deltav  =  dh_gmI(vz,0,vdah,vdah);

   deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl2,Nr2,Nt2,Nb2,socfd) -...
       (DRs_l*Ut2' + DRs_r*Ub2' + (Ul2')*DRs_t + (Ur2')*DRs_b)'); 
      deltav  =  Oh8_dh_gmI(vz,0,h,vdah,Bl2,Br2,Bt2,Bb2,socfd);
      for j=1:sx
          %M(1:sx+1:sx^2) = IDH(j)*e;
          M(1:sy+1:sy^2) = IDH(j)*e;
          VV(j,:) = (M\deltab(j,:)');
      end
        
    d2vz = (W*VV + deltav/h^2);
    if(ord == 2) d3vz = vdah; d4vz = vdah; d5vz = vdah; d6vz = vdah; return; end
    
    %==================  3rd 
    
    b =(2*al*bt*dvz.*vz + (bt-1)*dvz);
    
    deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl3,Nr3,Nt3,Nb3,socfd) - (DRs_l*Ut3' + DRs_r*Ub3' + (Ul3')*DRs_t + (Ur3')*DRs_b)'); 
        deltav  =  Oh8_dh_gmI(dvz,0,h,vdah,Bl3,Br3,Bt3,Bb3,socfd);

        for j=1:sx
            M(1:sy+1:sy^2) = IDH(j)*e;
            VV(j,:) = (M\deltab(j,:)');
        end

    d3vz = (W*VV + deltav/h^2);
    if(ord == 3) d4vz = vdah; d5vz = vdah;  d6vz = vdah; return; end
    
    %==================  4th  
    
    b =( 2*al*bt*(dvz.*dvz + vz.*d2vz) + (bt-1)*d2vz);
    
        deltab = W'*( Oh8_dh_gmI(b,0,h,vdah,Nl4,Nr4,Nt4,Nb4,socfd) - (DRs_l*Ut4' + DRs_r*Ub4' + (Ul4')*DRs_t + (Ur4')*DRs_b)'); 
        deltav  =  Oh8_dh_gmI(d2vz,0,h,vdah,Bl4,Br4,Bt4,Bb4,socfd);

        for j=1:sx
            M(1:sy+1:sy^2) = IDH(j)*e;
            VV(j,:) = (M\deltab(j,:)');
        end

    d4vz = (W*VV + deltav/h^2);
    if(ord == 4) d5vz = vdah;  d6vz = vdah; return; end
    
    %==================  5th
    
    b =( 2*al*bt*(3*dvz.*d2vz + vz.*d3vz) + (bt-1)*d3vz);

        deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl5,Nr5,Nt5,Nb5,socfd) - (DRs_l*Ut5' + DRs_r*Ub5' + (Ul5')*DRs_t + (Ur5')*DRs_b)'); 
        deltav  =  Oh8_dh_gmI(d3vz,0,h,vdah,Bl5,Br5,Bt5,Bb5,socfd);
        
        for j=1:sx
            M(1:sy+1:sy^2) = IDH(j)*e;
            VV(j,:) = (M\deltab(j,:)');
        end

    d5vz = (W*VV + deltav/h^2);
    if(ord == 5)  d6vz = vdah; return; end
        %==================  6th
    
    b =( 2*al*bt*(3*d2vz.*d2vz + 4*dvz.*d3vz + vz.*d4vz) + (bt-1)*d4vz);


        deltab = W'* (Oh8_dh_gmI(b,0,h,vdah,Nl6,Nr6,Nt6,Nb6,socfd) - (DRs_l*Ut6' + DRs_r*Ub6' + (Ul6')*DRs_t + (Ur6')*DRs_b)'); 
        deltav  =  Oh8_dh_gmI(d4vz,0,h,vdah,Bl6,Br6,Bt6,Bb6,socfd);

        for j=1:sx
            M(1:sy+1:sy^2) = IDH(j)*e;
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