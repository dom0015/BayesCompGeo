function [ He_n ] = hermite_poly_normalized( n )
%HERMITE_POLY Summary of this function goes here
%   Detailed explanation goes here
He_n=zeros(n);
He_n(1,1)=1;
He_n(2,2)=1;
diff=1:(n-1);
for i=3:n
    He_n(i,2:end)=He_n(i,2:end)+He_n(i-1,1:(end-1));
    He_n(i,1:(end-1))=He_n(i,1:(end-1))-diff.*He_n(i-1,2:end);
end
for i=1:n
    He_n(i,:)=He_n(i,:)/(sqrt(factorial(i-1)));
end
end

