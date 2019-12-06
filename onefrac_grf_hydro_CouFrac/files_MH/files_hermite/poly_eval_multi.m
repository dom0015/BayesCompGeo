function [ phi ] = poly_eval_multi( hermite, poly, x )
%POLY_EVAL Summary of this function goes here
%   Detailed explanation goes here
[N,dim]=size(poly);
phi=ones(N,1);
for i=1:N   % cyklus pres vsechny prvky polynomialni baze
    for j=1:dim     % cyklus pres jednotlive 1d polynomy
        % vsechny body v s1 pocitam najednou akorat po jednotlivych
        % dimenzich
        phi(i)=phi(i)*poly_eval(hermite(poly(i,j)+1,:),x(j));
    end
end