function [ lnratio_approx,Ap_approx,Bp,A_approx,Gp_approx ] = MCMC_alpha_approx( p, y, G_approx, noiseInvCov, priorInvCov, mu0, Au_approx, Bu )
%MCMC_ Summary of this function goes here
%   Detailed explanation goes here
%Am=norm(y-G(m))^2;
%Ap=norm(y-G(p))^2;
%Bm=norm(H(m)-mu0)^2;
%Bp=norm(p-mu0)^2;
Gp_approx=G_approx(p); % approximated observation operator in "p"
Ap_approx=(y-Gp_approx)'*noiseInvCov*(y-Gp_approx)/2;
Bp=(p-mu0)'*priorInvCov*(p-mu0)/2;
A_approx=Au_approx-Ap_approx;
lnratio_approx=A_approx+Bu-Bp;
end

