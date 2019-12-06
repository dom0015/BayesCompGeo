function [ values ] = poly_eval( p, grid )
%POLY_EVAL Summary of this function goes here
%   Detailed explanation goes here
n=find(abs(p)>0,1,'last');
values=zeros(size(grid));
temp=ones(size(grid));
values=values+p(1);
for i=2:n
    temp=temp.*grid;
    values=values+temp*p(i);
end
end