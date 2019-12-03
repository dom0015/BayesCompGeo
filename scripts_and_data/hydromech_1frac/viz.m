% load('uloha2.mat')
% load('hydromech_1frac.mat')
y=Q(2:end)';

no_windows=4;
sigma_num=1e-11;
sigma_obs=1e-13;
sigma_both=sigma_num^2/no_windows+sigma_obs^2;
noiseCov=diag(sigma_both*ones(no_windows,1));
noiseInvCov=inv(noiseCov);

likelihood=zeros(300,1);
for i=1:300
    Gu=ALL_Q(i,:)';
    likelihood(i)=-0.5*(y-Gu)'*noiseInvCov*(y-Gu);
end

prior=-0.5*(idx+6).^2;

figure;
plot(idx,exp(likelihood),'r')
hold on
plot(idx,exp(prior),'m')
plot(idx,exp(likelihood+prior),'b')

legend('likelihood', 'prior', 'posterior')