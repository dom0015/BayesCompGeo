function [ lnratio,Ap,Gp ] = MCMC_alpha_exact_dcg( y, Gp, noiseInvCov,Au,A_approx )
%MCMC_ Summary of this function goes here
%   Detailed explanation goes here
% Gp=G(p);
Ap=(y-Gp)'*noiseInvCov*(y-Gp)/2;
lnratio=Au-Ap-A_approx;
%disp(['lnratio=' num2str(lnratio)]);
end

