function [resmat,resmat0,TT,dy2]=deltah2d_Oh4(sx,sw)

if(nargin == 2)
    sw1 = 0;
end

    Mv=(-60/12)*diag(ones(sx,1));
    Mv(1,1)=-59/12;  Mv(1,2)=16/12;   Mv(1,3)=-1/12;
    Mv(2,1)=16/12;                    Mv(2,3)=16/12;  Mv(2,4)=-1/12;
    
    for l=3:sx-2
        Mv(l,l-2)=-1/12;
        Mv(l,l+2)=-1/12;
        Mv(l,l-1)=16/12;
        Mv(l,l+1)=16/12;
    end
          
    Mv(sx-1,sx-3)=-1/12; Mv(sx-1,sx-2)=16/12;                  Mv(sx-1,sx)=16/12;
                         Mv(sx,sx-2)=-1/12; Mv(sx,sx-1)=16/12;  Mv(sx,sx)=-59/12;
    
    
    resmat = Mv;
    resmat0 = Mv;
    resmat0(1,1)=-58/12;
    resmat0(sx,sx)=-58/12;
    for l=2:sx-1
        resmat0(l,l)=-59/12;
    end
  %{  
    
    TT = zeros(sx*sy);
    if(sw1==0)
        
        for i = 2:sy-1
            
        end
    end
%} 
    I=eye(sx);
    sy = sx;
    TT = zeros(sx*sy);
    TT(1:sx,1:sx) = resmat0;
    TT(sx*(sy-1)+1:end, sx*(sy-1)+1:end) = resmat0;
    for i=1:sy-2
       TT(i*sx+1:(i+1)*sx,i*sx+1:(i+1)*sx) = Mv;
    end
    
    for i=1:sy-2
       TT(i*sx+1:(i+1)*sx,i*sx+1-sx:(i+1)*sx - sx) = 4/3*I;
       TT(i*sx+1:(i+1)*sx,i*sx+1+sx:(i+1)*sx + sx) = 4/3*I;
    end
    TT((sy-1)*sx+1:(sy-1+1)*sx,(sy-1)*sx+1-sx:(sy-1+1)*sx - sx) = 4/3*I;
    TT(0*sx+1:(0+1)*sx,0*sx+1+sx:(0+1)*sx + sx) = 4/3*I;
    
    for i=2:sy-3
       TT(i*sx+1:(i+1)*sx,i*sx+1-2*sx:(i+1)*sx - 2*sx) = -1/12*I;
       TT(i*sx+1:(i+1)*sx,i*sx+1+2*sx:(i+1)*sx + 2*sx) = -1/12*I;
    end
    TT((sy-1)*sx+1:(sy-1+1)*sx,(sy-1)*sx+1-2*sx:(sy-1+1)*sx - 2*sx) = -1/12*I;
    TT(0*sx+1:(0+1)*sx,0*sx+1+2*sx:(0+1)*sx + 2*sx) = -1/12*I;
    
    dy2 = TT;
    dd = -2.5*ones(1,sx);
    dd(1) = -(2.5 - 1/12);
    dd(end) = dd(1);
    for i=0:sy-1
       dy2(i*sx+1:(i+1)*sx,i*sx+1:(i+1)*sx) = diag(dd);
    end
   