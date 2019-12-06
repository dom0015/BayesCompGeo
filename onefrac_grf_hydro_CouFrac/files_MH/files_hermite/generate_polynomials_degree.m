function [ poly ] = generate_polynomials_degree( dim, degree )
%GENERATE_POLYNOMIALS_DEGREE Summary of this function goes here
%   Detailed explanation goes here
if dim==1
    poly=(0:degree)';
   return; 
end
dim_idx=1:dim;
start=2;
poly=zeros(nchoosek(dim+degree,degree),dim);
for i=1:degree
%    temp= nchoosek(reshape(repmat(dim_idx,i,1),1,dim*i),i);
%    temp= unique(temp,'rows');
    temp=nmultichoosek(dim_idx,i);
   temp_poly=zeros(nchoosek(dim+i-1,i),dim);
   for j=1:dim
      temp_poly(:,j)=sum(temp==j,2);
   end
   
   
   poly(start:(start-1+nchoosek(dim+i-1,i)),:)= temp_poly;
   start=start+nchoosek(dim+i-1,i);
    
end


end

function combs = nmultichoosek(values, k)
%// Return number of multisubsets or actual multisubsets.
if numel(values)==1 
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
end