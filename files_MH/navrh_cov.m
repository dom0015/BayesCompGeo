function [cov_mat] = navrh_cov(sigma_veci,sigma_noise)
%NAVRH_COV Summary of this function goes here
%   Detailed explanation goes here

n=length(sigma_veci);
cov_mat_all=cell(n,1);

for i=1:n
cov_mat_all{i}=[sigma_veci(i)^2/4 -sigma_veci(i)^2/2+sigma_noise(i)^2/4;
    -sigma_veci(i)^2/2+sigma_noise(i)^2/4 sigma_veci(i)^2];
end
cov_mat=blkdiag(cov_mat_all{:});
cov_mat=cov_mat([(1:n)*2-1 (1:n)*2],[(1:n)*2-1 (1:n)*2]);
end
