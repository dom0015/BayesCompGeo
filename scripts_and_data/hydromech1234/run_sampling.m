% parallel.importProfile('/apps/all/MATLAB/2015b-EDU/SalomonPBSPro.settings')
% cluster = parcluster('local');
% pool=parpool(cluster,16);
%% Bayesian solution of an inverse problem
addpath(genpath('../../files_MH'))
addpath(genpath('../../files_elasticity'))
addpath(genpath('../../files_hydro'))
addpath(genpath('../../files_other'))

Nbatch=1;
no_windows=4;
chosen_confs = [1];
reduction=0.5;
initial_cov=[0.1,0.1];

maxtimeMH=0.2*60*10;%*60;
maxtimeDAMHSMU=0.4*60*10;%*60;
maxtimeDAMHSMU2=0.6*60*10;%*60;
maxtimeDAMH=0.95*60*10;%*60;

%% forward model and its parameters
no_confs=length(chosen_confs);
chosen_windows_str = num2str(sum(chosen_confs.*10.^(no_confs-1:-1:0)));

ulohy=cell(no_confs,1);
for i=chosen_confs
    ulohy{i}=load(['uloha' num2str(i) '.mat']);
end
G=@(x)G_operator(10.^x,ulohy);
% G = @(x)G_for_tests(x);

u_real=[-6 -6.5 -5.8 -7]';
Npar=length(u_real);

G_real=G(u_real);
Nobs=length(G_real);
disp(ceil(log10(abs(G_real))))

%% Bayes parameters
gamma=[1 1 1 1]'; % prior standard deviation
priorInvCov=diag(1./(gamma.^2));
mu0=[-6 -6 -6 -6]';     % prior mean
sigma_num=1e-11;
sigma_obs=1e-13;
sigma_both=sigma_num^2/no_windows+sigma_obs^2;
% sigma_both=1;
noiseCov=diag(sigma_both*ones(no_windows,1)); temp2=noiseCov;
noiseInvCov=inv(noiseCov); temp=noiseInvCov;
for i=2:no_confs
    noiseCov=blkdiag(noiseCov,temp2);
    noiseInvCov=blkdiag(noiseInvCov,temp);
end

rng(2)
noise=chol(noiseCov)'*randn(Nobs,1);
y=G_real+noise;
save(['noisy_observation_' num2str(no_windows) 'w.mat'],'y','noise','G_real');

% load(['noisy_observation_' num2str(no_windows) 'w.mat'])
% y=G_real;
% pattern=false(4*no_windows,1);
% for i=1:no_confs
%     pattern(((chosen_confs(i)-1)*no_windows+1):chosen_confs(i)*no_windows)=1;
% end
% y=y(pattern);

%% sampling process parameters
propCov=navrh_cov([1 1],initial_cov);
% propCov=eye(4);
propMH=chol(propCov*reduction^2)';
% propDAMHSMU=propMH;
% propDAMHSMU2=chol(propCov*(reduction*1.5)^2)';
% propDAMH=propDAMHSMU2;
stepsMH=10000;
stepsDAMH=1000000;
stepsDAMHSMU=1000000;
stepsDAMHSMU2=1000000;
period_savingMH=100;
period_savingDAMHSMU=200;
period_savingDAMHSMU2=200;
period_savingDAMH=200;
% sigmaMH=0.005*gamma;    % standard deviation of the proposal density
% sigmaDAMHSMU=0.005*gamma;    % standard deviation of the proposal density
% sigmaDAMHSMU2=0.03*gamma;    % standard deviation of the proposal density
% sigmaDAMH=0.03*gamma;    % standard deviation of the proposal density

%% surrogate model parameters
approx_type=1;      % approximation type - 0 rbf, 1 col
SCM_maxdegree=7;

label=[num2str(no_windows) 'w_' num2str(Npar) '_' num2str(length(y)) '_' chosen_windows_str ];
%% surrogate model update
if approx_type==0
    kernel=@(x)x.^3;
    upgrade_approx=@(s1,sG)upgrade_approx_rbf( s1, sG, kernel);
    label = [label '_rbf'];
else
    upgrade_approx=@(s1,sG,sM)upgrade_approx_col_compress( s1, sG, sM, Npar, SCM_maxdegree);
    label = [label '_col'];
end

time_start=tic;
[status, msg, msgID] = mkdir('res');
save(['res/CONFIGURATION0_' label '.mat'],'-v7.3');

LAST_SAMPLE=cell(Nbatch,1);
ALL_U_cell=cell(Nbatch,1);
ALL_GU_cell=cell(Nbatch,1);
ALL_S_cell=cell(Nbatch,1);
ALL_M_cell=cell(Nbatch,1);
recent_cov_cell=cell(Nbatch,1);

for batch=1:Nbatch
    label_par=[label '_' num2str(batch)];
    rng(batch+100); u0_MH=my_normrnd(mu0,gamma/1000);
    u0_MH=u_real;
    [SAMPLES_MH,OBSERVATIONS_MH,MULTIPLICITY_MH,EVALUATED_U,EVALUATED_GU, TIMES, W,~,~,recent_cov ]...
        = MCMC_standard_dcg_compress(u0_MH,y,G,noiseInvCov,priorInvCov,mu0,propMH,stepsMH,...
        period_savingMH,time_start,maxtimeMH,batch+200,label_par,initial_cov,reduction);
    LAST_SAMPLE{batch}=SAMPLES_MH(end,:)';
    ALL_U_cell{batch}=EVALUATED_U; %%
    ALL_GU_cell{batch}=EVALUATED_GU; %%
    ALL_S_cell{batch}=SAMPLES_MH;
    ALL_M_cell{batch}=MULTIPLICITY_MH;
    recent_cov_cell{batch}=recent_cov;
end

ALL_U=cell2mat(ALL_U_cell);
ALL_GU=cell2mat(ALL_GU_cell);
NU=size(ALL_U,1);
ALL_S=cell2mat(ALL_S_cell);
ALL_M=cell2mat(ALL_M_cell);
recent_cov=mean(cell2mat(recent_cov_cell),1);
% propCov=my_cov(ALL_S,ALL_M);
% propDAMHSMU=chol(propCov*1)';
propDAMHSMU = update_proposal_cov(recent_cov,ALL_S,ALL_M,reduction);
propDAMHSMU2 = update_proposal_cov(recent_cov,ALL_S,ALL_M,reduction*1.5);
propDAMH=propDAMHSMU2;

% my_corr_plot(ALL_S,ALL_M,u_real,prior_mean,prior_std)
% cov_generate(navrh_cov([1 1],recent_cov))
for batch=1:Nbatch
    label_par=[label '_' num2str(batch)];
    [SAMPLES0,OBSERVATIONS0,MULTIPLICITY0, ~, log_approx, EVALUATED_U,EVALUATED_GU ]...
        = MCMC_approx_dcg_compress(LAST_SAMPLE{batch},y,G,noiseInvCov,priorInvCov,mu0,propDAMHSMU,stepsDAMHSMU,...
        period_savingDAMHSMU,time_start,maxtimeDAMHSMU,batch+300,[label_par '_prep'],...
        ALL_U,ALL_GU,ones(NU,1),upgrade_approx,approx_type );
    LAST_SAMPLE{batch}=SAMPLES0(end,:)';
    ALL_U_cell{batch}=EVALUATED_U;
    ALL_GU_cell{batch}=EVALUATED_GU;
end

ALL_U=[ALL_U; cell2mat(ALL_U_cell)];
ALL_GU=[ALL_GU; cell2mat(ALL_GU_cell)];
NU=size(ALL_U,1);

for batch=1:Nbatch
    label_par=[label '_' num2str(batch)];
    [SAMPLES0,OBSERVATIONS0,MULTIPLICITY0, G_approx, log_approx, EVALUATED_U,EVALUATED_GU ]...
        = MCMC_approx_dcg_compress(LAST_SAMPLE{batch},y,G,noiseInvCov,priorInvCov,mu0,propDAMHSMU2,stepsDAMHSMU2,...
        period_savingDAMHSMU2,time_start,maxtimeDAMHSMU2,batch+300,label_par,...
        ALL_U,ALL_GU,ones(NU,1),upgrade_approx,approx_type );
    LAST_SAMPLE{batch}=SAMPLES0(end,:)';
    ALL_U_cell{batch}=EVALUATED_U;
    ALL_GU_cell{batch}=EVALUATED_GU;
end

ALL_U=[ALL_U; cell2mat(ALL_U_cell)];
ALL_GU=[ALL_GU; cell2mat(ALL_GU_cell)];
NU=size(ALL_U,1);
G_approx=upgrade_approx(ALL_U,ALL_GU,ones(NU,1));

for batch=1:Nbatch
    label_par=[label '_' num2str(batch)];
    [SAMPLES,OBSERVATIONS,MULTIPLICITY, G_approx_rbf, log_approx, EVALUATED_U,EVALUATED_GU ]...
        = DAMH_dcg(LAST_SAMPLE{batch},y,G,noiseInvCov,priorInvCov,mu0,propDAMH,stepsDAMH,...
        period_savingDAMH,time_start,maxtimeDAMH,batch+400,label_par,G_approx );
    LAST_SAMPLE{batch}=SAMPLES(end,:)';
    ALL_U_cell{batch}=EVALUATED_U;
    ALL_GU_cell{batch}=EVALUATED_GU;
end

ALL_U=[ALL_U; cell2mat(ALL_U_cell)];
ALL_GU=[ALL_GU; cell2mat(ALL_GU_cell)];
NU=size(ALL_U,1);
G_approx=upgrade_approx(ALL_U,ALL_GU,ones(NU,1));

save(['res/CONFIGURATION1_' label '.mat'],'-v7.3');

% viz_corrplot(Nbatch,label,u_real,gamma,mu0);
% load(['MH_ALL_' label '.mat'])
% my_corr_plot(MH_SAMPLES,MH_MULTIPLICITY,u_real,prior_mean,prior_std)
% figure; corrplot(ALL_U)
% 
% function [ Gx, other ] = G_operator(x)
%     Gx=[x(1)*x(2); sin(3*x(3)); x(4)];
%     other=0;
% end
