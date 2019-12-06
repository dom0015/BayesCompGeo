%% Bayesian solution of an inverse problem
addpath(genpath('files_MH'))

%% necessary parameters
n_windows=2; % number of windows on the outflow side
chosen_windows = [1 3]; n_configs = length(chosen_windows);
pattern = repmat([false true(1,n_windows)],1,n_configs);
s=2;                % number of random parameters (length of u)
batch=14;
gamma=2;            % prior standard deviation
mu0=[-7 -7]';       % prior mean
sigma=1e-11*ones((n_windows+1)*n_configs,1);
sigma(pattern)=sigma(pattern)/n_windows;
u0=mu0;             % initial iteration of the MH algorithm
steps_MH = 200;
steps_MH_approx=1e6;         % number of all MH steps
period_saving=1e4;
sigmaMH=gamma*ones(s,1);    % standard deviation of the instrumental density
approx_type=1;      % approximation type - 0 rbf, 1 col

%% observation operator G
[ A,b,freeNode,U0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac ] = preparation_windows(n_windows);
G=@(x)observation_windows( A,b,freeNode,U0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,chosen_windows,x,n_windows );

%% information from data y = G(k) + gaussian noise
u_real=[-6 -8.5]';
G_real=G(u_real);
rng(3)
noise=my_normrnd(zeros(length(G_real),1),sigma);
y=G_real+noise;

%% rbf approximation of G
if approx_type==0
    kernel=@(x)x.^3;
    upgrade_approx=@(s1,sG)upgrade_approx_rbf( s1, sG, kernel);
else
    upgrade_approx=@(s1,sG)upgrade_approx_col( s1, sG, s, 5);
end

%% Metropolis-Hastings with approximation
time_start=tic;
label=['two13' '_' num2str(s) '_' num2str(length(G_real)) '_' num2str(steps_MH_approx)];
if approx_type==0
    label = [label '_rbf'];
else
    label = [label '_col'];
end
save(['res/CONFIGURATION0_' label '.mat'],'-v7.3');
%visualization( mu0, gamma, u_real ); hold on
[evaluated_u_col,evaluated_Gu_col,evaluated_u_rbf,evaluated_Gu_rbf]=MCMC_standard( u0,y,G,sigma,gamma,mu0,0.1*sigmaMH,steps_MH,100,batch*2,label);
time_MH=toc(time_start);
if approx_type==0
    evaluated_u_init = evaluated_u_rbf;
    evaluated_Gu_init = evaluated_Gu_rbf;
else
    evaluated_u_init = evaluated_u_col;
    evaluated_Gu_init = evaluated_Gu_col;
end
visualization( mu0, gamma, u_real ); hold on
[G_approx_rbf,log_approx,evaluated_Gu_mh]=MCMC_approx(evaluated_u_col(end,:)',y,G,sigma,gamma,mu0,sigmaMH,steps_MH_approx,period_saving,batch,label,evaluated_u_init,evaluated_Gu_init,upgrade_approx,approx_type);
time_MH_approx=toc(time_start)-time_MH;
save(['res/CONFIGURATION_' label '.mat'],'-v7.3');




