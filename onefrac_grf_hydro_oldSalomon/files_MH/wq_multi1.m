function [ w,q ] = wq_multi1( kernel, s1, Gapp, G )
%WQ_MULTI1 Summary of this function goes here
%   Detailed explanation goes here
% RBF for changes in only the last coordinate,
% but with multiple "last coordinates".
% "vstupem" jsou souradnice [a b x1], [a b x2], [a b x3]
%              a souradnice [a b y1], [a b y2], [a b y3]
% vystupem jsou vektory         wq1       wq2       wq3

[n_s,m]=size(s1); % pocet ridicich bodu, pocet parametru
d=m+1; % pocet parametru plus jedna
n_G=size(Gapp,2);
w=zeros(n_s,n_G);
q=zeros(d+1,n_G);
S1=zeros(n_s);
for i=1:m
    temp=repmat(s1(:,i),1,n_s);
    S1=S1+(temp-temp').^2;
end

for i=1:n_G
    s=[s1 Gapp(:,i)];
    s_=[s1 G(:,i)];
    temp=repmat(Gapp(:,i),1,n_s);
    SG=(temp-temp').^2;
    S=kernel(sqrt(S1+SG));
    P=[s ones(n_s,1)];
    A=[S P; P' zeros(d+1)];
    B=[s_-s; zeros(d+1,d)];

    wq=gmres(A,B(:,d),100,1e-6,100);
    w(:,i)=wq(1:n_s,:);
    q(:,i)=wq(n_s+1:end,:);
end

end

