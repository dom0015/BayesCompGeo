function [ G_approx ] = upgrade_approx_col_compress( s1, sG, sM, dim, max_degree)%, y, sigma, gamma, mu0 )
%UPGRADE_APPROX_COL Summary of this function goes here
%   dim ... pocet nahodnych parametru
Q=size(s1,1);   % pocet kolokacnich bodu
%dim3=[4    10    20    35    56    84   120   165   220   286];
%degree=ceil(log(Q/(dim*dim)))+1;
degree=floor(log(Q)/log(dim));
if degree==0
    degree=1;
end
disp(degree)
%degree=floor(log(Q))
if max_degree<degree
    degree=max_degree;
end
% degree=max_degree;
poly=generate_polynomials_degree(dim, degree);
N=size(poly,1);     % velikost polynomialni baze
hermite=hermite_poly_normalized(degree+1); % tabulka koeficientu vsech 1d hermitovskych polynomu do urciteho radu
% do poly_eval dosazujeme radky hermite(i,:)

coefs=ones(Q,N);
for i=1:N   % cyklus pres vsechny prvky polynomialni baze
    for j=1:dim     % cyklus pres jednotlive 1d polynomy
        % vsechny body v s1 pocitam najednou akorat po jednotlivych
        % dimenzich
        coefs(:,i)=coefs(:,i).*poly_eval(hermite(poly(i,j)+1,:),s1(:,j));
    end
end
coefs_M=coefs.*repmat(sM,1,N);
A=coefs_M'*coefs;
%disp([length(s1) N rank(coefs) rank(A) ])
b=coefs_M'*sG;
c = pinv(A)*b;
% c=A\b;
% c=gmres(A,b,1000,1e-8,1000);
G_approx=@(x)c'*poly_eval_multi(hermite,poly,x);
end

