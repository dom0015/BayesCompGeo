%% Bayesian solution of an inverse problem
addpath(genpath('files_MH'))

%% necessary parameters
n_windows=4; % number of windows on the outflow side
chosen_windows = [1]; n_configs = length(chosen_windows);
pattern = repmat([false true(1,n_windows)],1,n_configs);
s1=4;                % number of random parameters (length of u)
s=1*(s1+0+1);
batch=14;
gamma=1;            % prior standard deviation
mu0=0*ones(s,1);       % prior mean
mu0(1:2)=1e-6;
sigma=1e-11*ones((n_windows+1)*n_configs,1);
sigma(pattern)=sigma(pattern)/n_windows;
u0=mu0;             % initial iteration of the MH algorithm
stepsMH = 200;
stepsDAMH=10000;         % number of all MH steps
period_savingMH=100;
period_savingDAMH=1000;
sigmaMH=0.5*gamma*ones(s,1);    % standard deviation of the instrumental density
sigmaMH(1)=1e-8;
sigmaDAMH=0.5*gamma*ones(s,1);    % standard deviation of the instrumental density
sigmaDAMH(1)=1e-8;
approx_type=1;      % approximation type - 0 rbf, 1 col
SCM_maxdegree=5;

%% observation operator G
[A,b,freeNode,U0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,fracture_positions ] = preparation_windows(n_windows);
no_frac=length(fracture_positions);
for i=1:no_frac
    [data_generator{i}, ~] = set_fracture(s1,fracture_positions{i},intersections(:,i) );
end
G=@(x)observation_windows( A,b,freeNode,U0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,x(1),chosen_windows,x(2:end),data_generator,n_windows,0 );
G_visual=@(x)observation_windows( A,b,freeNode,U0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,x(1),chosen_windows,x(2:end),data_generator,n_windows,1 );

%% information from data y = G(k) + gaussian noise
u_real=[1e-6];
for i=1:no_frac
    u_real=[u_real; data_generator{i}.random_coef];
end
G_real=G(u_real);
%G_visual(u_real);
rng(3)
noise=my_normrnd(zeros(length(G_real),1),sigma);
y=G_real;%+noise;

% rbf approximation of G
if approx_type==0
    kernel=@(x)x.^3;
    upgrade_approx=@(s1,sG)upgrade_approx_rbf( s1, sG, kernel);
else
    upgrade_approx=@(s1,sG,sM)upgrade_approx_col_compress( s1, sG, sM, s, SCM_maxdegree);
    upgrade_approx2=@(s1,sG,sM)upgrade_approx_col( s1, sG, s, SCM_maxdegree);
end

%% Metropolis-Hastings with approximation
time_start=tic;
label=['twofrac_grf' '_' num2str(s) '_' num2str(length(G_real)) '_' num2str(stepsDAMH)];
if approx_type==0
    label = [label '_rbf'];
else
    label = [label '_col'];
end
save(['res/CONFIGURATION0_' label '.mat'],'-v7.3');
%visualization( mu0, gamma, u_real ); hold on

% [evaluated_u_init2,evaluated_Gu_init2,evaluated_u_rbf,evaluated_Gu_rbf]=                                           MCMC_standard( u0,y,G,sigma,gamma,mu0,sigmaMH,stepsMH,100,batch*2,label);
[evaluated_u_col,evaluated_Gu_col,evaluated_M_init,evaluated_u_rbf,evaluated_Gu_rbf, TIMES, W ] = MCMC_standard_dcg_compress( u0,y,G,sigma,gamma,mu0,sigmaMH,stepsMH,period_savingMH,batch*2,label);
time_MH=toc(time_start);
if approx_type==0
    evaluated_u_init = evaluated_u_rbf;
    evaluated_Gu_init = evaluated_Gu_rbf;
else
    evaluated_u_init = evaluated_u_col;
    evaluated_Gu_init = evaluated_Gu_col;
end
% visualization( mu0, gamma, u_real ); hold on
% [G_approx_rbf,log_approx,evaluated_Gu_mh]=                                                     MCMC_approx(evaluated_u_col(end,:)',y,G,sigma,gamma,mu0,sigmaDAMH,stepsDAMH,period_savingDAMH,batch,label,evaluated_u_init2,evaluated_Gu_init2,upgrade_approx2,approx_type);

[SAMPLES,OBSERVATIONS,MULTIPLICITY, G_approx_rbf, log_approx, evaluated_Gu_mh ] = MCMC_approx_dcg_compress(evaluated_u_col(end,:)',y,G,sigma,gamma,mu0,sigmaDAMH,stepsDAMH,period_savingDAMH,batch,label,evaluated_u_init(2:end,:),evaluated_Gu_init(2:end,:),evaluated_M_init(2:end),upgrade_approx,approx_type );


time_MH_approx=toc(time_start)-time_MH;
save(['res/CONFIGURATION_' label '.mat'],'-v7.3');

load(['res/MH_', label, '.mat'])
figure; plot(SAMPLES','.b'); hold on; plot(u_real,'*r')
