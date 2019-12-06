function [ vG ] = apply_multi1( v1,kernel,w,q,s1,Gapp )
%APPLY_MULTI1 Summary of this function goes here
%   Detailed explanation goes here

%% application of the resulting function to all nodes
n=size(v1,1);        % no of nodes to be morphed
[n_s,m]=size(s1); % pocet ridicich bodu, pocet parametru
n_G=size(Gapp,2);
N1=zeros(n_s,n);   
for i=1:m
    temp=s1(:,i)*ones(1,n)-ones(n_s,1)*v1(:,i)';
    temp=temp.^2;
    N1=N1+temp;
end
vG=zeros(n,n_G);
for i=1:n_G
    %temp=Gapp(:,i)*ones(1,n)-ones(n_s,1)*vGapp(:,i)';
    %NG=temp.^2;
    %N=sqrt(N1+NG);              % distances between all points in "s" and in "v"
    N = sqrt(N1);
    temp1=kernel(N)'*w(:,i);     % first sum from d
    %temp2=[v1 vGapp(:,i) ones(n,1)]*q(:,i);	% second sum from d
    temp2=[v1 zeros(n,1) ones(n,1)]*q(:,i);	% second sum from d
    displacement=temp1+temp2;
    vG(:,i)=vG(:,i)+displacement;
    
end

