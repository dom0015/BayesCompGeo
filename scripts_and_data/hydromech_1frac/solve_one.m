
addpath(genpath('../../files_MH'))
addpath(genpath('../../files_elasticity'))
addpath(genpath('../../files_hydro'))
addpath(genpath('../../files_other'))

%% PROBLEM ASSEMBLY
Nxy=101;
L1=10; L2=L1;
sumbdomains_FETI=100;

mat_const=1e9;
frac_press_val=1;
frac_start_end={[0.1 0.5], [0.9 0.5]};

mat_omega_const=1e-15;
mat_frac_const=1e-6;
alfa_inter_const=1e-5;
cislo_ulohy=2;
n_windows=4;
         
addpath(genpath('files_elasticity'))
addpath(genpath('files_hydro'))
FETI_problem_assembly
tocouple_aperture_independent

SMALSE_params.rel=1.0e-9;
SMALSE_params.rho0=1;
SMALSE_params.betarho=2;
SMALSE_params.Gama = 1;
SMALSE_params.M_start=0.5;
SMALSE_params.tol_to_update=1e3;
SMALSE_params.maxiter_cg = 2000;
SMALSE_params.type='m';
SMALSE_params.print=true;
SMALSE_params.print_couple=true;
SMALSE_params.coupling_iter=100; 
SMALSE_params.eps_coupling=1e-6;
problem_setting.B_i=[];
problem_setting.lambda_ker=[];

clearvars -except SMALSE_params problem_setting hydro_problem
save(['uloha_test.mat'],'problem_setting','hydro_problem','SMALSE_params')

% load('uloha2.mat')
% SMALSE_params.print=true;
% SMALSE_params.print_couple=true;

param=10.^(-4.5);

%% COUPLED SOLVER
hydro_problem.mat_frac=param;
% depends on fracture aperture:
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
d = 1e-4*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
D_old=D;

res_D=[];
res_PRESSURE=[];
alpha=0.5;
% beta=1e-1;
rel_D_=[];
rel_P_=[];
for i=1:SMALSE_params.coupling_iter
    
    % HYDRO
    [PRESSURE,~,ugrad,Q,PRESSURE_diff]=tocouple_handle(D,hydro_problem);
    res_PRESSURE(:,i)=cell2mat(PRESSURE);
    
    if i>1
        for j=1:size(PRESSURE_diff,1)
            for jj=1:size(PRESSURE_diff,2)
                PRESSURE_diff{j,jj}=(PRESSURE_diff{j,jj}*alpha+(1-alpha)*PRESSURE_old{j,jj});
            end
        end
        ugrad=ugrad*alpha+(1-alpha)*ugrad_old;
    end
%     for j=1:size(PRESSURE_diff,1)
%         for jj=1:size(PRESSURE_diff,2)
%             PRESSURE_diff{j,jj}=(PRESSURE_diff{j,jj}*beta);
%         end
%     end
%     ugrad=beta*ugrad;
    
    % MECH
    [problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE_diff,-1*ugrad);
    [D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params);
    
    for j=1:length(D)
        D{j}=max(D{j},1e-10);
    end
    
    res_D(:,i)=cell2mat(D);
    if i>1
        for j=1:length(D)
            D{j}=(D{j}*alpha+(1-alpha)*D_old{j});
        end
    end
    
    tmp=sqrt(sum((res_PRESSURE(:,1:end-1)-res_PRESSURE(:,2:end)).^2,1))./sqrt(sum((res_PRESSURE(:,2:end)).^2,1));
    tmp2=sqrt(sum((res_D(:,1:end-1)-res_D(:,2:end)).^2,1))./sqrt(sum((res_D(:,2:end)).^2,1));
    
    D_old=D;
    PRESSURE_old=PRESSURE_diff;
    ugrad_old=ugrad;
    
    % SMALSE precision to match coupled precision
%     if i>2
%         SMALSE_params.rel=tmp2(end)/1000;
%     end
    % coupling parameters update
%     alpha=min(0.8,alpha^(0.5));
%     disp(["alpha", num2str(alpha)]);
%     beta=min(1,beta^(0.5)+1e-1);
    if i>1
        if alpha>=0.1 && tmp(end)<SMALSE_params.eps_coupling && tmp2(end)<SMALSE_params.eps_coupling
            break
        end
    end
    
    if i>1
        rel_D_=[rel_D_ tmp(end)];
        rel_P_=[rel_P_ tmp2(end)];
    end
end

figure(5)
subplot(1,2,1)
hold on
plot(rel_P_,'k-','LineWidth',1)
set(gca,'YScale','log')
%figure(5)
subplot(1,2,2)
yyaxis right
hold on
plot(rel_D_,'k-','LineWidth',1)
set(gca,'YScale','log')

if i>1
    fprintf('i=%d: alpha = %.2f || pres = %d || d = %d\n',i,alpha,tmp(end),tmp2(end));
end
