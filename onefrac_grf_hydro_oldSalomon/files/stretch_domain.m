function [ node ] = stretch_domain( node )
%STRETCH_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

idx=node(:,1)>9;
node(idx,1)=(node(idx,1)-9).^6*1000+9;

idx=node(:,2)>9;
node(idx,2)=(node(idx,2)-9).^6*1000+9;

idx=node(:,1)<1;
node(idx,1)=-(1-node(idx,1)).^6*1000+1;

idx=node(:,2)<1;
node(idx,2)=-(1-node(idx,2)).^6*1000+1;

end

