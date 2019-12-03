% parallel.importProfile('/apps/all/MATLAB/2015b-EDU/SalomonPBSPro.settings')
% cluster = parcluster('local');
% pool=parpool(cluster,24);
%% Bayesian solution of an inverse problem
addpath(genpath('files_MH'))
addpath(genpath('files_elasticity'))
addpath(genpath('files_hydro'))
addpath(genpath('files_other'))
% load('SCM.mat');

% par
for batch=1
    tic
    %% forward model and its parameters
    ulohy=cell(4,1);
    for i=1:4
        ulohy{i}=load(['uloha' num2str(i) '.mat']);
    end
    %params=[1 1.5 1 1]*10^(-6);
    G=@(x)G_operator(10.^x,ulohy);
    u_real=[-6 -6.5 -5.8 -7]';
    [G_real,aperture]=G(u_real);
    Npar=4;
    Nobs=16;
    toc

    %% sampling process parameters
    gamma=1;                % prior standard deviation
    mu0=-6*ones(Npar,1);     % prior mean
    sigma=0.01*abs(G_real);   % noise std
    sigma=0.01*mean(abs(G_real))*ones(Nobs,1);   % noise std
    sigma=sqrt(((1e-11)^2)/4+(1e-13)^2)*ones(Nobs,1);
    stepsMH = 100;
    stepsDAMH=1000;
    period_savingMH=10;
    period_savingDAMH=100;
    sigmaMH=0.05*gamma*ones(Npar,1);    % standard deviation of the proposal density
    sigmaDAMH=0.1*gamma*ones(Npar,1);    % standard deviation of the proposal density

    %% surrogate model parameters
    approx_type=1;      % approximation type - 0 rbf, 1 col
    SCM_maxdegree=7;

    %G_approx=@(x)c'*poly_eval_multi(hermite,poly,x);

    %% information from data y = G(k) + gaussian noise
    rng(3)
    noise=my_normrnd(zeros(length(G_real),1),sigma);
    y=G_real+noise;
    
    %% surrogate model update
    if approx_type==0
        kernel=@(x)x.^3;
        upgrade_approx=@(s1,sG)upgrade_approx_rbf( s1, sG, kernel);
    else
        upgrade_approx=@(s1,sG,sM)upgrade_approx_col_compress( s1, sG, sM, Npar, SCM_maxdegree);
    end 

    %% Metropolis-Hastings with approximation
    time_start=tic;
    label=['twofrac_grf' '_' num2str(Npar) '_' num2str(length(G_real)) '_' num2str(stepsDAMH) '_' num2str(batch)];
    if approx_type==0
        label = [label '_rbf'];
    else
        label = [label '_col'];
    end
    save(['res/CONFIGURATION0_' label '.mat'],'-v7.3');
    rng(batch+200)
    [evaluated_u_col,evaluated_Gu_col,evaluated_M_init,evaluated_u_rbf,evaluated_Gu_rbf, TIMES, W ] = MCMC_standard_dcg_compress( normrnd(mu0,gamma),y,G,sigma,gamma,mu0,sigmaMH,stepsMH,period_savingMH,batch+100,label);
    timeMH=toc(time_start);
    if approx_type==0
        evaluated_u_init = evaluated_u_rbf;
        evaluated_Gu_init = evaluated_Gu_rbf;
    else
        evaluated_u_init = evaluated_u_rbf;
        evaluated_Gu_init = evaluated_Gu_rbf;
    end
    rng(batch+300)
    [SAMPLES,OBSERVATIONS,MULTIPLICITY, G_approx_rbf, log_approx, evaluated_Gu_mh ] = MCMC_approx_dcg_compress(evaluated_u_col(end,:)',y,G,sigma,gamma,mu0,sigmaDAMH,stepsDAMH,period_savingDAMH,batch,label,evaluated_u_init(2:end,:),evaluated_Gu_init(2:end,:),evaluated_M_init(2:end),upgrade_approx,approx_type ); 
    %[SAMPLES,OBSERVATIONS,MULTIPLICITY, G_approx_rbf, log_approx, evaluated_Gu_mh ] = DAMH_dcg(u_real,y,G,sigma,gamma,mu0,sigmaDAMH,stepsDAMH,period_savingDAMH,batch+400,label,G_approx );

%     timeDAMH=toc(time_start)-timeMH;
%     save(['res/CONFIGURATION1_' label '.mat'],'-v7.3');
end

% load(['res/MH_', label, '.mat'])
% figure; plot(SAMPLES','.b'); hold on; plot(u_real,'*r')
% load(['res/DAMH_', label, '.mat'])
% plot(SAMPLES','.g');
