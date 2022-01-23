function [UUp,PUp,thetaVector,solutionNorms,tauVector,angl,sw_div]=... 
    sol_ch_1d_v6(U,x,t,prmtrs,bt1,bt2,al,c,theta,derivative,P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    sw = 0;  
    if (nargin == 11) 
        sw=1; 
    end

    % constants
    tau = prmtrs.tau;
    h=prmtrs.h;
    eps = prmtrs.eps;
    tauMax = 10*tau;
    tauIncreasedIteration = 60;
    tauDecreasedIteration = 1;
    iterMax = prmtrs.iterMax;
    bt = bt1/bt2;
    autoStop = 0;
    afterCounter = 5000;
    [zeroX,zeroT]=GetZeroNodes(x,t);
    
    manualStop=0;
    boundaryHit = 0;
    sw_div = 0;

    UvsUupInfNorm = ones(1,iterMax);
    thetaVector = zeros(1,iterMax);
    tauVector = zeros(1,iterMax); 
    tauVector(1) = tau;
    residualInfNorm = ones(1,iterMax/10);
    angl = zeros(1,iterMax/10);
    step = Step(h);
    zeroMatrix = zeros(size(U));
    
    Utt = YDerivative(U, zeroMatrix, derivative.second);
    Uxx = XDerivative(U, zeroMatrix, derivative.second);
    if(sw == 0)
        P = theta*al*bt*U.^2 - Uxx + bt * Utt;        
    end
    
    sx = length(x);
    st = length(t);
    dF1 = zeros( size(U) ); 
    dF2 = zeros( size(U) ); 
    for jt=2:st-1
        for ix = 2:sx-1
            dU1l = -bt/h^2;
            dU1c = (-2*bt/tau^2 + 2*bt/h^2 );
            dU1r = -bt/h^2;
            dU1d = bt/tau^2;
            dU1u = bt/tau^2;
            dP1l = -1/h^2;
            dP1c = 2/h^2;
            dP1r = -1/h^2;
            dF1(ix,jt) = dU1l * U(ix-1,jt) + dU1c * U(ix,jt) + dU1r * U(ix+1,jt) +...
                dU1d * U(ix,jt-1) + dU1u * U(ix,jt+1) + dP1l * P(ix-1,jt) + dP1c * P(ix,jt) + dP1r * P(ix+1,jt);
            
            dU2l = -1/h^2;
            dU2c = (-2*bt/tau^2 + 2/h^2 + 2*al*bt*U(ix,jt) );
            dU2r = -1/h^2;
            dU2d = dU1d;
            dU2u = dU1u;
            dP2c = -1;
            dF2(ix,jt) = dU2l * U(ix-1,jt) + dU2c * U(ix,jt) + dU2r * U(ix+1,jt) +...
                dU2d * U(ix,jt-1) + dU2u * U(ix,jt+1) + dP2c * P(ix,jt);
        end
    end
end

