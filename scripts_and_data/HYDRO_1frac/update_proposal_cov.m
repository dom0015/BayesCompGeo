function [propL,recent_cov] = update_proposal_cov(recent_cov,a,b,reduction)
%UPDATE_PROPOSAL_COV Summary of this function goes here
%   Detailed explanation goes here
propCov=eye(11);
propL=chol(propCov*reduction^2)';
end

