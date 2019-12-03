function [ lnratio,Ap,Bp,Gp ] = MCMC_alpha_dcg( p, y, Gp, noiseInvCov, priorInvCov, mu0, Au, Bu )
%MCMC_ Summary of this function goes here
%   Detailed explanation goes here
%Am=norm(y-G(m))^2;
%Ap=norm(y-G(p))^2;
%Bm=norm(H(m)-mu0)^2;
%Bp=norm(p-mu0)^2;
% Gp=G(p); % approximated observation operator in "p"
Ap=(y-Gp)'*noiseInvCov*(y-Gp)/2;
Bp=(p-mu0)'*priorInvCov*(p-mu0)/2;
lnratio=Au-Ap+Bu-Bp;
% disp(Au)
% disp(-Ap)
% disp(Bu)
% disp(-Bp)
% disp(lnratio)
end

