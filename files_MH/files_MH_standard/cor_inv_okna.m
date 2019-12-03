function [out] = cor_inv_okna(k,sigma,sigma2)
%COR_INV_OKNA Summary of this function goes here
%   Detailed explanation goes here
M=sigma^2*ones(k)/k^2+eye(k)*sigma2^2;
out=inv(M);
% [U,D,V]=svd(M);
% out=U(:,1)*D(1)^(-1)*V(:,1)';
end

