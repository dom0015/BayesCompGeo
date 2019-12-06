function [ G_approx_rbf ] = upgrade_approx_rbf( s1, sG, kernel )
%UPGRADE_APPROX_RBF Summary of this function goes here
%   Detailed explanation goes here
[w,q]=wq_multi1( kernel, s1, 0*sG, sG );
G_approx_rbf=@(x)apply_multi1(x',kernel,w,q,s1,0*sG)';

end

