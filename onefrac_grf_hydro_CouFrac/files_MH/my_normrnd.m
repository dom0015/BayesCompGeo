function [ r ] = my_normrnd( mu, sigma )
%MY_NORMRND Summary of this function goes here
%   Detailed explanation goes here
[mu_rows,mu_cols]=size(mu);
[sigma_rows,sigma_cols]=size(sigma);
if mu_rows<sigma_rows
    ratio=sigma_rows/mu_rows;
    mu=repmat(mu,ratio,1);
elseif mu_rows>sigma_rows
    ratio=mu_rows/sigma_rows;
    sigma=repmat(sigma,ratio,1); 
end
if mu_cols<sigma_cols
    ratio=sigma_cols/mu_cols;
    mu=repmat(mu,1,ratio);
elseif mu_cols>sigma_cols
    ratio=mu_cols/sigma_cols;
    sigma=repmat(sigma,1,ratio); 
end
r=randn(size(mu,1),size(mu,2)).*sigma+mu;

end

