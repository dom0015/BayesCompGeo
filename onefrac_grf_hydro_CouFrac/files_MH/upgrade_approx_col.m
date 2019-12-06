function [ G_approx ] = upgrade_approx_col( s1, sG, dim, degree)%, y, sigma, gamma, mu0 )
%UPGRADE_APPROX_COL Summary of this function goes here
%   dim ... pocet nahodných parametru
%   POUZE PRO JEDNU DIMENZI G:Rn->R1
poly=generate_polynomials_degree(dim, degree);
N=size(poly,1);     % velikost polynomialni baze
hermite=hermite_poly_normalized(degree+1); % tabulka koeficientu vsech 1d hermitovskych polynomu do urciteho radu
% do poly_eval dosazujeme radky hermite(i,:)
Q=size(s1,1);   % pocet kolokacnich bodu
coefs=ones(Q,N);
for i=1:N   % cyklus pres vsechny prvky polynomialni baze
    for j=1:dim     % cyklus pres jednotlive 1d polynomy
        % vsechny body v s1 pocitam najednou akorat po jednotlivych
        % dimenzich
        coefs(:,i)=coefs(:,i).*poly_eval(hermite(poly(i,j)+1,:),s1(:,j));
    end
end
% zbytecne... pro vyhodnoceni pi(s1) nepotrebujeme cyklus
%post=exp(-(y-sG).^2/(2*sigma*sigma)-sqrt(sum((s1-mu0).*(s1-mu0),2))/(2*gamma*gamma))

% A=zeros(N);
% b=zeros(N,1);
% for i=1:N
%     b(i)=dot(coefs(:,i),sG);
%     for j=1:i
%         if i==j
%             A(i,i)=dot(coefs(:,i),coefs(:,i));
%         else
%             A(j,i)=dot(coefs(:,j),coefs(:,i));
%             A(i,j)=A(j,i);
%         end
%     end
% end
A=coefs'*coefs;
% disp([length(s1) rank(A) rank(coefs)])
b=coefs'*sG;
% c=A\b;
c = pinv(A)*b;
% c=gmres(A,b,1000,1e-8,1000);
G_approx=@(x)c'*poly_eval_multi(hermite,poly,x);
end

