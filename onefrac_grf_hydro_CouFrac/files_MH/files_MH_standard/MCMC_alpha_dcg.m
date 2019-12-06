function [ lnratio,Ap,Bp,Gp ] = MCMC_alpha_dcg( p, y, Gp, sigma, gamma, mu0, Au, Bu )
%MCMC_ Summary of this function goes here
%   Detailed explanation goes here
%Am=norm(y-G(m))^2;
%Ap=norm(y-G(p))^2;
%Bm=norm(H(m)-mu0)^2;
%Bp=norm(p-mu0)^2;
% Gp=G(p); % approximated observation operator in "p"
Ap=(y-Gp)'*diag(1./(sigma.^2))*(y-Gp)/2;
Bp=(p-mu0)'*diag(1./(gamma.^2))*(p-mu0)/2;
lnratio=Au-Ap+Bu-Bp;
end

