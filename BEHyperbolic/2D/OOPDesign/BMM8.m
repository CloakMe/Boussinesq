function [Y]=BMM8(ma3x,X)
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
    
    if(~(sm == 9))
        error('ERROR; size of >ma3x< should be 9!,i.e. a cross section of the band! ');
    end
    sx = size(X,2);
    ml = size(X,1);
    
    if(sm == 3)
        Y(1) = ma3x(2:3)*X(1:2)';
        %for l=2:sx-1
            Yc=ma3x*[X(1:end-2); X(2:end-1); X(3:end)];
        %end
        Ysx = ma3x(1:2)*X(sx-1:sx)';
        Y = [Y(1) Yc Ysx];
    else
        Y(1) = ma3x(5:9)*X(1:5)';
        Y(2) = ma3x(4:9)*X(1:6)';
        Y(3) = ma3x(3:9)*X(1:7)';
        Y(4) = ma3x(2:9)*X(1:8)';
        %for l=3:sx-2
            %Y(l)=ma3x*X(l-2:l+2);
            Yc=ma3x*[X(1:end-8); X(2:end-7); X(3:end-6); X(4:end-5); X(5:end-4);
                                                     X(6:end-3); X(7:end-2); X(8:end-1); X(9:end)];
        %end
        Ysx3 = ma3x(1:8)*X(sx-7:sx)';
        Ysx2 = ma3x(1:7)*X(sx-6:sx)';
        Ysx1 = ma3x(1:6)*X(sx-5:sx)';
        Ysx =  ma3x(1:5)*X(sx-4:sx)';
        Y = [Y(1) Y(2) Y(3) Y(4) Yc Ysx3 Ysx2 Ysx1 Ysx];
    end
end
