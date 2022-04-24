function [Y]=BMM(ma3x,X,m11)
%Function for band matrix multiplication
%F(M) = ma3x; Y = M*X; m11 := M(1,1).
% >ma3x< and >m11< contain the band matrix structure.
%M is 5 diagonal or 3 diagonal symmetric matrix where each band of M
%is represented by only one real number i.e. band = a*[1 1 ... 1]; 
%Exception for 5 diag matrix M: The main diagonal (band) equals: c * [(m11/c) 1 1 ... 1 (m11/c)];
%F(M) is the crossection of the band, i.e. 
% [a b a] if M is 3 diag or
% [a b c b a ] if M is 5 diag.

sm = size(ma3x,1)*size(ma3x,2);
if(nargin == 2)
    if(sm == 3)
        m11 = -345; % if sm == 3 this value does not matter
    else
        error('ERROR; missmatch of >ma3x<`size or M(1,1) not given! ');
    end
end
    
    if(sm == 5 || sm == 3)
    else
        error('ERROR; size of >ma3x< should be 5 or 3!,i.e. a cross section of the band! ');
    end
    sx = size(X,2);
    if(sx<3)
        error('ERROR; X should be a ROW vector or size of X too small! ');
    end
    
    if(sm == 3)
        Y(1) = ma3x(2:3)*X(1:2)';
        %for l=2:sx-1
            Yc=ma3x*[X(1:end-2); X(2:end-1); X(3:end)];
        %end
        Ysx = ma3x(1:2)*X(sx-1:sx)';
        Y = [Y(1) Yc Ysx];
    else
        Y(1) = m11*X(1) + ma3x(4:5)*X(2:3)';
        Y(2) = ma3x(2:5)*X(1:4)';
        %for l=3:sx-2
            %Y(l)=ma3x*X(l-2:l+2);
            Yc=ma3x*[X(1:end-4); X(2:end-3); X(3:end-2); X(4:end-1); X(5:end)];
        %end
        Ysx1 = ma3x(1:4)*X(sx-3:sx)';
        Ysx =  ma3x(1:2)*X(sx-2:sx-1)' + m11*X(sx);
        Y = [Y(1) Y(2) Yc Ysx1 Ysx];
    end
