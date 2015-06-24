classdef BEUtilities
  methods(Static)
      
    function [c,T,errOrder]=GetFinDiffCoeff(nodesPosition, derivativeOrder)

        % Determine the coefficents to a finite differencing scheme
        % ie. find the coeffs [c1,c2,c3] that are appropriate for the 2nd order
        %     cent. diff. approx. to the first derivative :

        % (du/dx)[j] = 1/dx*(c1*u[j-1] + c2*u[j] + c3*u[j+1]) + error terms

        % Inputs:
        %     x - vector containing relative distances to sampling points
        %       ie. for example above would be x=[-1,0,1]
        %       ie. 5 terms full forward: x = [0,1,2,3,4]
        %       ie. 5 terms central :     x = [-2,-1,0,1,2]
        %
        %     derv - order of the derivative to approximate

        % Ouputs:
        %     c - the coefficients you're looking for!
        %     T optional output matrix showing the Taylor Table coefficents
        %     errOrder - lowest order of higher order terms
        %                (order of accuracy, ie. H.O.T. = O(dx^2) is 2nd order) 
        %The code:
        N=numel(nodesPosition);  %number of unknowns
        T=zeros(N,N);
        for j=1:N  % do for each x
            for n=0:N+6 % Do for each coeff of the Taylor expansion
                        % (Expand out 6 extra terms for higher order accuracies
                        %   at higher derivatives)
                T(j,n+1)=-1/factorial(n)*nodesPosition(j)^n;      
            end
        end
        % Vector holding place of the derivative
        d=zeros(N,1);
        d(derivativeOrder+1)=-1;

        % Calculate the coefficients
        c=T(:,1:N)'\d;

        for i=1:numel(c)  % set any small value to 0
            if abs(c(i))<1E-10
                c(i)=0;
            end
        end

        % Determine the order of accuracy


        a=c'*T; %  Multiply the Taylor table coefficients by the finite difference
                %    coefficients and add up all like terms
                %  Terms that don't add to zero are the objective derivative
                %    and the residual higher order terms

            for i = 1:numel(a) % make sure values near zero = zero
                if abs(a(i))<1E-5;
                    a(i)=0;
                end
            end

        pos=find(a,2);    % location of first nonzero term is the order of
                          % derivative you're approximating. Second nonzero
                          % term tells the lowest order of the higher order terms
        errOrder=pos(2)-1-derivativeOrder;
                          % ie.   H.O.T. = O(dx^2) is 2nd order accuracy
                          %       H.O.T. = O(dx^4) is 4th order accuracy

        T=T(:,1:N); % Only output Taylor Table terms that are used  
    end
    
    function x=PentSolv(a11,A,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        % Solve a pentadiagonal system Dx=b where D is a strongly nonsingular
        % matrix, D is symmetric and each subdiagonal has the following vector form:
        % a*ones(1,N-l) == [a a .... a], (l=0,1,2  depending on the diagonal).
        % The main diagonal has the following form: [a11 b b .... b a11]
        % Here the input matrix ( A ) presents the matrix D in a shortened way
        % A = D(3,1:5)! :
        %      1 2 3 4 5 6 
        %    -------------------------------------
        % 1  | # $ &                             |
        % 2  | $ * $ &                           |
        % 3  | & $ * $ &                         | = D
        % 4  |   & $ * $ &                       |
        % 5  |           ...                     |
        %
        % i.e.: A = [ & $ * $ & ] and a11 = #
        %
        % If D is not a pentadiagonal matrix, results will be wrong
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(size(A,1)*size(A,2)~=5)
            error('A should be a vector, of type [a b c b a]!');
        end
        [M,N]=size(b);
        if((M>1 && N==1) || (M==1 && N>1))
            if(M>N)
                N=M;
            end
        else
            error('b should be a vector, not number or matrix!');
        end

        x=zeros(N,1);

            % Extract bands
            d=A(3);
            f=A(2);
            e=A(1);

            alpha=zeros(N,1);
            gamma=zeros(N-1,1);
            delta=zeros(N-2,1);
            z=zeros(N,1);

            % Factor A=LDL'
            alpha(1)=a11;
            gamma(1)=f/alpha(1);
            delta(1)=e/alpha(1);

            alpha(2)=d-f*gamma(1);
            gamma(2)=(f-e*gamma(1))/alpha(2);
            delta(2)=e/alpha(2);

            for k=3:N-2
                alpha(k)=d-e*delta(k-2)-alpha(k-1)*gamma(k-1)^2;
                gamma(k)=(f-e*gamma(k-1))/alpha(k);
                delta(k)=e/alpha(k);
            end

            alpha(N-1)=d-e*delta(N-3)-alpha(N-2)*gamma(N-2)^2;
            gamma(N-1)=(f-e*gamma(N-2))/alpha(N-1);
            alpha(N)=a11-e*delta(N-2)-alpha(N-1)*gamma(N-1)^2;

            % Update Lx=b, Dc=z

            z(1)=b(1);
            z(2)=b(2)-gamma(1)*z(1);

            for k=3:N
                z(k)=b(k)-gamma(k-1)*z(k-1)-delta(k-2)*z(k-2);
            end

            c=z./alpha;

            % Backsubstitution L'x=c
            x(N)=c(N);
            x(N-1)=c(N-1)-gamma(N-1)*x(N);

            for k=2:N-1
                x(N-k)=c(N-k)-gamma(N-k)*x(N-k+1)-delta(N-k)*x(N-k+2);
            end
    end

    function x=TridiagSolv(A,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Solve a tridiagonal system Dx=b where D is a strongly nonsingular
        % matrix, D is symmetric and each subdiagonal has the following vector form:
        % a*ones(1,N-l) == [a a .... a], (l=0,1,2  depending on the diagonal).
        % The main diagonal has the following form: [a11 b b .... b a11]
        % Here the input matrix ( A ) presents the matrix D in a shortened way
        % A = D(2,1:3)! :
        %      1 2 3 4 5 6 
        %    ------------------------------------
        % 1  | # $                              |
        % 2  | $ * $                            |
        % 3  |   $ * $                          | = D
        % 4  |     $ * $                        |
        % 5  |          ...                     |
        %
        % i.e.: A = [ $ * $ ] and a11 = #
        % If D is not a 3 diagonal matrix, results will be wrong!
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(size(A,1)*size(A,2)~=3)
            error('A should be a vector, of type [b c b]!');
        end
        [M,N]=size(b);
        if((M>1 && N==1) || (M==1 && N>1))
            if(M>N)
                N=M;
            end
        else
            error('b should be a vector, not number or matrix!');
        end


        x=zeros(N,1);

        % Extract bands
        d=A(2);
        f=A(1);

        alpha=zeros(N,1);
        gamma=zeros(N-1,1);
        z=zeros(N,1);

        % Factor A=LDL'
        alpha(1)=d;
        gamma(1)=f/alpha(1);

        alpha(2)=d-f*gamma(1);
        gamma(2)=(f)/alpha(2);

        for k=3:N-2
            alpha(k)=d-alpha(k-1)*gamma(k-1)^2;
            gamma(k)=(f)/alpha(k);
        end

        alpha(N-1)=d-alpha(N-2)*gamma(N-2)^2;
        gamma(N-1)=(f)/alpha(N-1);
        alpha(N)=d-alpha(N-1)*gamma(N-1)^2;

        % Update Lx=b, Dc=z

        z(1)=b(1);
        z(2)=b(2)-gamma(1)*z(1);

        for k=3:N
            z(k)=b(k)-gamma(k-1)*z(k-1);
        end

        c=z./alpha;

        % Backsubstitution L'x=c
        x(N)=c(N);
        x(N-1)=c(N-1)-gamma(N-1)*x(N);

        for k=2:N-1
            x(N-k)=c(N-k)-gamma(N-k)*x(N-k+1);
        end
 
    end  
    
    function [Y]=BandMatMult(ma3x,X,m11)
        %Function for band matrix multiplication
        %F(M) = ma3x; Y = M*X; m11 := M(1,1).
        % >ma3x< and >m11< contain the band matrix structure.
        %M is 5 diagonal or 3 diagonal symmetric matrix where each band of M
        %is represented by only one real number i.e. band = a*[1 1 ... 1]; 
        %Exception for 5 diag matrix M: The main diagonal (band) equals: c * [(m11/c) 1 1 ... 1 (m11/c)];
        %F(M) is the crossection of the band, i.e. 
        % [a b a] if M is 3 diag or
        % [a b c b a ] if M is 5 diag.
        
        if(nargin == 2)
            Y = BEUtilities.BandMatMult2(ma3x,X);
            return;
        end
        
        sm = size(ma3x,1)*size(ma3x,2);
        if(sm == 5 || sm == 3)
        else
            error('ERROR; size of >ma3x< should be 5 or 3!,i.e. a cross section of the band! ');
        end
        sx = size(X,2);
        ml = size(X,1);

        if(ml>1 && sx~=1)
           if(ml<5)
               size(X)
             error('BMM; right side is MATRIX; MATRIX size too small small (<5)!!! ');
           end
           if(sm == 3)
          % HM = zeros(ml-2,sx);
           HM = X(1:end-2,:);
           XV(1,:) = HM(:)';
           HM = X(2:end-1,:);
           XV(2,:) = HM(:)';
           HM = X(3:end,:);
           XV(3,:) = HM(:)';   

           Y = zeros(ml,sx);

             Y(1,:) = ma3x(2:3)*X(1:2,:);

               Yc=vec2mat(XV'*ma3x',ml-2)';

             Ysx = ma3x(1:2)*X(end-1:end,:);

            Y = [Y(1,:); Yc; Ysx];
           else
                error('ERROR; size(ma3x) must equal [1,3] ... not ready yet for 5 diag matrices! ');
           end
            return;
        end

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
    end
    
    function [dx2,O4dx2]=GetFinDiffMat(sx,h)
        %
        % Input:
        % sx - size of the matrix 
        % h - discretization step size
        % Output:
        % First matrix is second order finite diff matrix
        % Second matrix is fourth order finite diff matrix
        %               %
        dx2=-2*diag(ones(sx,1));
        dx2(1,2)=1;
        for l=2:sx-1
            dx2(l,l-1)=1;
            dx2(l,l+1)=1;
        end
        dx2(sx,sx-1)=1;
        %dx2 = dx2/h^2;

        O4dx2=-2.5*diag(ones(sx,1));
        O4dx2(1,1)=-2.5;   O4dx2(1,2)=16/12; O4dx2(1,3)=-1/12;
        O4dx2(2,1)=16/12;                    O4dx2(2,3)=16/12; O4dx2(2,4)=-1/12;
        for l=3:sx-2
            O4dx2(l,l-2)=-1/12;
            O4dx2(l,l+2)=-1/12;
            O4dx2(l,l-1)=16/12;
            O4dx2(l,l+1)=16/12;
        end
        O4dx2(sx-1,sx-3)=-1/12; O4dx2(sx-1,sx-2)=16/12;                       O4dx2(sx-1,sx)=16/12;    
                                O4dx2(sx,sx-2)=-1/12; O4dx2(sx,sx-1)=16/12;   O4dx2(sx,sx)=-2.5;
        return;
    end
    
    function [sdah] = GetDerOnBnd(X,Y,sdah,t,bt,c,mu,ord)

        c12 = 1 - c^2;
        cbt=c*sqrt(bt);

        sdah(:,:,1) = -6 * mu * (c12 - 1 + cbt ^ 2) * (c12^2 * X.^4 + (-6 * c12 *(cbt*t - Y).^2) .* X.^2 +...
            (cbt * t - Y).^4) ./ (c12 * X.^2 + (Y - cbt * t).^2).^4;
        if(ord == 2)
            return;
        end

        sdah(:,:,2) = 24 * mu * (cbt .* t - Y) .* cbt .* (cbt ^ 2 + c12 - 1) .* (5 * c12 ^ 2 * X.^4 -...
            10 .* c12 .* (cbt * t - Y).^2 .* X.^2 + (cbt * t - Y).^4) ./ (c12 * X.^2 + (cbt * t - Y).^2).^ 5;
        if(ord == 3)
            return;
        end

        sdah(:,:,3) = -120.* mu.* cbt ^ 2.* (cbt ^ 2 + c12 - 1).* (-c12.* X.^2 + (cbt.* t - Y).^2) .* (c12 ^ 2.* X.^4 -...
            14.* c12.* (cbt.* t - Y).^2.* X.^2 + (cbt.* t - Y).^4) ./ (c12.* X.^2 + (cbt.* t - Y).^2).^ 6;
        if(ord == 4)
            return;
        end

        sdah(:,:,4) = 720.* cbt ^ 3.* mu.* (cbt.* t - Y).* (cbt ^ 2 + c12 - 1) .* (-7.* c12 ^ 3.* X.^6 +...
            35.* c12 ^ 2.* (cbt.* t - Y).^2.* X.^4 - 21.* c12.* (cbt.* t - Y).^4.* X.^2 +...
            (cbt.* t - Y).^6) ./ (c12.* X.^2 + (cbt.* t - Y).^2).^ 7;
        if(ord == 5)
            return;
        end

        sdah(:,:,5) = -5040.* mu.* cbt ^ 4.* (cbt ^ 2 + c12 - 1).* (c12 ^ 4.* X.^8 - 28.* c12 ^ 3.* (cbt.* t - Y).^2.* X.^6 +...
            70.* c12 ^ 2.* (cbt.* t - Y).^4.* X.^4 - 28.* c12.* (cbt.* t - Y).^6.* X.^2 +...
            (cbt.* t - Y).^8) ./ (c12.* X.^2 + (cbt.* t - Y).^2).^ 8;
        if(ord == 6)
            return;
        end

        sdah(:,:,6) = 40320.* (cbt.* t - Y).* mu.* cbt ^ 5.* (cbt ^ 2 + c12 - 1).* (-3.* c12.* X.^2 +...
            (cbt.* t - Y).^2).* (-3.* c12 ^ 3.* X.^6 + 27.* c12 ^ 2.* (cbt.* t - Y).^2.* X.^4 -...
            33.* c12.* (cbt.* t - Y).^4.* X.^2 + (cbt.* t - Y).^6) ./ (c12.* X.^2 + (cbt.* t - Y).^2).^ 9;

        if(ord == 7)
            return;
        end

        sdah(:,:,7) = -362880.* cbt ^ 6.* mu.* (cbt ^ 2 + c12 - 1).* (-c12.* X.^2 + (cbt.* t - Y).^2).* (c12 ^ 4.* X.^8 -...
            44.* c12 ^ 3.* (cbt.* t - Y).^2.* X.^6 + 166.* c12 ^ 2.* (cbt.* t - Y).^4.* X.^4 -...
            44.* c12.* (cbt.* t - Y).^6.* X.^2 + (cbt.* t - Y).^8) ./ (c12.* X.^2 + (cbt.* t - Y).^2).^ 10;

    end

  end
  
  methods (Access = private, Static)
      
    function [Y]=BandMatMult2(ma3x,X)
        %Function for band matrix multiplication
        %F(M) = ma3x; Y = M*X; 
        % >ma3x< and >m11< contain the band matrix structure.
        %M is 5 diagonal or 3 diagonal symmetric matrix where each band of M
        %is represented by only one real number i.e. band = a*[1 1 ... 1]; 
        %F(M) is the crossection of the band, i.e. 
        % ma3x = [a b a] if M is 3 diag or
        % ma3x = [a b c b a ] if M is 5 diag.

        sm = size(ma3x,1)*size(ma3x,2);
    
        if( sm ~= 5 && sm ~= 3 )
            error('ERROR; size of >ma3x< should be 5 or 3!,i.e. a cross section of the band! ');
        end
        sx = size(X,2);
        sy = size(X,1);

        if(sy>1 && sx~=1)
           if(sy<5)
               size(X)
             error('BMM; right side is MATRIX; MATRIX size too small small (<5)!!! ');
           end
           if(sm == 3)
               HM = X(1:end-2,:);
               XV(1,:) = HM(:)';
               HM = X(2:end-1,:);
               XV(2,:) = HM(:)';
               HM = X(3:end,:);
               XV(3,:) = HM(:)';   

               Y=vec2mat(XV'*ma3x',sy-2)';

           else
               HM = X(1:end-4,:);
               XV(1,:) = HM(:)';
               HM = X(2:end-3,:);
               XV(2,:) = HM(:)';
               HM = X(3:end-2,:);
               XV(3,:) = HM(:)';
               HM = X(4:end-1,:);
               XV(4,:) = HM(:)';
               HM = X(5:end-0,:);
               XV(5,:) = HM(:)';  

               Y=vec2mat(XV'*ma3x',sy-4)';
           end
            return;
        end
        
        if(sm == 3)
            if(sx<3)
                error('X must be a ROW vector! or size of X too small! ');
            end
            Y=ma3x*[X(1:end-2); X(2:end-1); X(3:end)];
        else
            if(sx<5)
                error('X must be a ROW vector! or size of X too small! ');
            end
            Y=ma3x*[X(1:end-4); X(2:end-3); X(3:end-2); X(4:end-1); X(5:end)];
        end
    end
    
  end
end