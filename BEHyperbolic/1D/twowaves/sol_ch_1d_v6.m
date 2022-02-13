function [UUp,PUp,thetaVector,solutionNorms,tauVector,angl,sw_div]=... 
    sol_ch_1d_v6(U,x,t,prmtrs,bt1,bt2,al,c,theta,finiteDiff,P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    sw = 0;  
    if (nargin == 11) 
        sw=1; 
    end
    
    derivatived2T = finiteDiff.secondT;
    derivatived2X = finiteDiff.secondX;
    % constants
    btExt = prmtrs.btExt;
    tau = prmtrs.tau;
    h=prmtrs.h;
    tTau = prmtrs.tTau;
    eps = prmtrs.eps;
    tauMax = 10*tau;
    tauIncreasedIteration = 60;
    tauDecreasedIteration = 1;
    iterMax = prmtrs.iterMax;
    bt = bt1/bt2;
    autoStop = 0;
    afterCounter = 5000;
    %[zeroX,zeroT]=GetZeroNodes(x,t);
    
    manualStop=0;
    boundaryHit = 0;
    sw_div = 0;

    UvsUupInfNorm = ones(1,iterMax);
    thetaVector = zeros(1,iterMax);
    tauVector = zeros(1,iterMax); 
    tauVector(1) = tau;
    residualInfNorm = ones(1,iterMax/10);
    angl = zeros(1,iterMax/10);
    %step = Step(h);
    zeroMatrix = zeros(size(U));
    
    Utt = YDerivativeEvenFunctions2D(U, zeroMatrix, derivatived2T);
    Uxx = XDerivativeEvenFunctions2D(U, zeroMatrix, derivatived2X);
    if(sw == 0)
        P = theta*al*bt*U.^2 - Uxx + bt * Utt; %      
    end
    
    sx = length(x);
    st = length(t);
    dF1 = zeros( size(U) ); 
    dF2 = zeros( size(U) );
    changeU = ( (U) );
    changeP = ( (P) );
    
    dU1l = -bt/h^2;
    dU1c = (-2*bt/tTau^2 + 2*bt/h^2 );
    dU1r = -bt/h^2;
    dU1d = bt/tTau^2;
    dU1u = bt/tTau^2;

    dP1l = -1/h^2;
    dP1c = 2/h^2;
    dP1r = -1/h^2;

    dU2l = -1/h^2;
    dU2r = -1/h^2;
    dU2d = dU1d*btExt;
    dU2u = dU1u*btExt;

    dP2c = -1;
    tic
    for it=2:st-1
        for jx = 2:sx-1
            dF1(jx,it) = dU1l * changeU(jx-1,it) + dU1c * changeU(jx,it) + dU1r * changeU(jx+1,it) +...
                dU1d * changeU(jx,it-1) + dU1u * changeU(jx,it+1) + dP1l * changeP(jx-1,it) + dP1c * changeP(jx,it) + dP1r * changeP(jx+1,it);
            
            dU2c = ( -2*bt*btExt/tTau^2 + 2/h^2 + theta*al*bt*U(jx,it) );    %        
            dF2(jx,it) = dU2l * changeU(jx-1,it) + dU2c * changeU(jx,it) + dU2r * changeU(jx+1,it) +...
                 dU2d * changeU(jx,it-1) + dU2u * changeU(jx,it+1) + dP2c * changeP(jx,it); % 
        end
    end  
    toc
      
    L=-GenerateBandDiagMatrix(derivatived2X, sx);
    M=bt*GenerateBandDiagMatrix(derivatived2T, st);
    
    P = btExt * U * M' + L * U + theta*al*bt*U.^2; %
    
    tic
    FF1 =          U * M' +  L * (bt * U + P);
    FF2 =  btExt * U * M' +  L * U - P + theta*al*bt*U .^ 2; %  
    toc
    
%     figure(5);
%     mesh( x(3:end-2), t(3:end-2), FF1(3:end-2,3:end-2)' );
%     xlabel('x');    ylabel('t');
%     title('FF1');
    
    tic
    dG1 =         bt * YDerivativeEvenFunctions2D(U, zeroMatrix, finiteDiff.secondT) - XDerivativeEvenFunctions2D( bt * U + P, zeroMatrix, finiteDiff.secondX);
    dG2 = btExt * bt * YDerivativeEvenFunctions2D(U, zeroMatrix, finiteDiff.secondT) - XDerivativeEvenFunctions2D(U, zeroMatrix, finiteDiff.secondX) + ...
        2*theta*al*bt*U .* U - P;
    toc
    [X,D] = eig(M);
    Lstar = zeros( sx, sx, st );
    subU = U;
    
    Isx = diag( ones(1,sx) );
    b = ( - FF1 ) * X;
    tic
    for it=1:st
        Ln = 2*theta*al*bt*diag(U(:,it));
        Lstar(:,:,it) = L * ( ( bt - btExt * D(it,it) ) * Isx - L - Ln) + D(it,it) *Isx ; %
        subU(:,it) = Lstar(:,:,it)\b(:,it);
    end
    shiftU = subU/X;
    toc
    shiftP = - FF2 - L * shiftU - btExt * shiftU * M' - 2*theta*al*bt*U .* shiftU; % 
    figure(4);
    mesh(x(2:end-1),t(2:end-1),shiftU(2:end-1,2:end-1)');
    xlabel('x');    ylabel('t');
    title('shift U');
%     figure(5);
%     mesh(x(2:end-1),t(3:end-2),shiftP(2:end-1,3:end-2)');
% 	  xlabel('x');    ylabel('t');
% 	  title('shift P');
%     figure(1)
%     mesh(x,t,dF1');
%     xlabel('x');    ylabel('t');
%     figure(2)
%     mesh(x(2:end-1),t(2:end-1),(dF1(2:end-1,2:end-1)-dG1(2:end-1,2:end-1))');
%     xlabel('x');    ylabel('t');
%     title('dF1-dG1');
%     figure(3)
%     mesh(x(3:end-2),t(2:end-1),(dF1(3:end-2,2:end-1)-dFF1(3:end-2,2:end-1))');
%     xlabel('x');    ylabel('t');
%     title('dF1-dFF1');
%     figure(3)
%     mesh(x,t,(dFF1-dF1)');
%     xlabel('x');    ylabel('t');
     
%     figure(5)
%     mesh(x(2:end-1),t(2:end-1),(dF2(2:end-1,2:end-1)- dG2(2:end-1,2:end-1))');
%     xlabel('x');    ylabel('t');
%     title('dF2-dG2');
%     figure(6)
%     mesh(x(2:end-1),t(2:end-1),(dF2(2:end-1,2:end-1)- dFF2(2:end-1,2:end-1))');
%     xlabel('x');    ylabel('t');
%     title('dF2-dFF2');
%     figure(6)
%     mesh(x,t,(dF2- dFF2)');
%     xlabel('x');    ylabel('t');
     UUp = U + shiftU;
     PUp = P + shiftP;
     thetaVector = theta;
     solutionNorms = 1;
     tauVector =1; 
     angl = 0;
     sw_div = 0;
end

