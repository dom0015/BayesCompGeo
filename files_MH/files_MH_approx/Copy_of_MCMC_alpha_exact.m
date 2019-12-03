function [ lnratio,Ap,Gp ] = MCMC_alpha_exact( p, y, G, sigma,Au,A_approx )
%MCMC_ Summary of this function goes here
%   Detailed explanation goes here
Gp=G(p);
Ap=(y-Gp)'*diag(1./(sigma.^2))*(y-Gp)/2;
lnratio=Au-Ap-A_approx;
%disp(['lnratio=' num2str(lnratio)]);
end

