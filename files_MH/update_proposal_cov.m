function [propMH,new_cov] = update_proposal_cov(recent_cov,SAMPLES,MULTIPLICITY,reduction)
%UPDATE_PROPOSAL_COV Summary of this function goes here
%   Detailed explanation goes here
w=0.5;
sigma2=1;
[~,COR,~]=my_cov(SAMPLES,MULTIPLICITY);
n13=max(COR(1,3),-1);
n24=max(COR(2,4),-1);
n13=sqrt((n13+1)*sigma2*2);
n24=sqrt((n24+1)*sigma2*2);
n13=n13*w+recent_cov(1)*(1-w);
n24=n24*w+recent_cov(2)*(1-w);

n13=max(n13,0.002);
n24=max(n24,0.002);
n13=min(n13,sqrt(2));
n24=min(n24,sqrt(2));
disp([n13 n24]);
new_cov=[n13 n24];
propCov = navrh_cov([1 1],new_cov);
propMH=chol(propCov*reduction^2)';
end

