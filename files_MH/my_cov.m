function [covs,cors,means] = my_cov(values,weights)
%MY_COV Summary of this function goes here
%   Detailed explanation goes here
[~,n]=size(values);

means=my_mean(values,weights);
covs=zeros(n);
for i=1:n
    for j=1:i
        covs(i,j)=my_Exy(values,weights,i,j)-means(i)*means(j);
        covs(j,i)=covs(i,j);
    end
end
cors=covs;
for i=1:n
    for j=1:n
        cors(i,j)=cors(i,j)/sqrt(covs(i,i))/sqrt(covs(j,j));
    end
end

end

function [m_w]=my_mean(values,weights)
n=size(values,2);
m_w=sum(values.*repmat(weights,1,n),1)/sum(weights);
end

function [xy_w]=my_Exy(values,weights,i,j)
xy_w=sum(values(:,i).*values(:,j).*weights)/sum(weights);
end
