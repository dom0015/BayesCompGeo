%% Bayesian solution of an inverse problem

%% necessary parameters
s=3;                % number of random parameters (length of u)
batch=14;
gamma=1;            % prior standard deviation
mu0=[8 6 4]';       % prior mean
sigma=10;           % noise standard deviation
u0=mu0;             % initial iteration of the MH algorithm
steps_MH = 100;
steps_MH_approx=1000;         % number of all MH steps
period_saving=10;
sigmaMH=1*ones(s,1);    % standard deviation of the instrumental density
approx_type=0;      % approximation type - 0 rbf, 1 col

%% observation operator G
G=@(x)[x(1)*x(1)*x(2); x(2)*x(3); x(2)*x(3); x(2)*x(3)];

%% information from data y = G(k) + gaussian noise
u_real=[9 5 3]';
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
label=['nic' '_' num2str(s) '_' num2str(length(G_real)) '_' num2str(steps_MH_approx)];
if approx_type==0
    label = [label '_rbf'];
else
    label = [label '_col'];
end
save(['res/CONFIGURATION_' label '.mat']);
visualization( mu0, gamma, u_real ); hold on
[evaluated_u_col,evaluated_Gu_col,evaluated_u_rbf,evaluated_Gu_rbf]=MCMC_standard( u0,y,G,sigma,gamma,mu0,0.1,steps_MH,period_saving,batch*2,label);
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
save(['res/CONFIGURATION_' label '.mat']);
