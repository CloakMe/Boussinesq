function [qD] = transf2qD(TF,x,y,ij,lo)
sx = length(x);
sy = length(y);
   RD = zeros(length(lo:sy-1),length(1:ij));
   RU = zeros(length(lo:sy),length(1:ij));
   RU = TF';
   RD = RU(2:end,:);
   RD = flipud(RD);

   R = [RD; RU];
   clear('RD');
   clear('RU');
%figure(3);   mesh(x(ij:sx),y,R(end-sy+1:end,:));
   L = zeros(size(R,1),size(R,2)-1);
   L = R(:,2:end);
   L = fliplr(L);
   %u_t0 = [L R]'; 
   %qD = u_t0(end-sy+1:end,:)';
   qD = [L(end-sy+1:end,:) R(end-sy+1:end,:)]'; 
  
   %{
   figure(12)
   mesh(x,y,qD');
   %mesh(x,y,dudt_t0');
   title('IC')
   xlabel('x');            ylabel('y');
   %}