
addpath(genpath('files_elasticity'))
addpath(genpath('files_hydro'))
%% load everything
ulohy=cell(4,1);
for i=1:4
    ulohy{i}=load(['uloha' num2str(i) '.mat']);
end
t=tic;
for ress=1:1

for c_ulohy=1:4
    problem_setting=ulohy{c_ulohy}.problem_setting;
    hydro_problem=ulohy{c_ulohy}.hydro_problem;
    SMALSE_params=ulohy{c_ulohy}.SMALSE_params;
    
    
    SMALSE_params.eps_coupling=1e-4;
    params=[1 1 1 1]*10^(-6);
    
    
    [Q,D] = coupled_solver(params,hydro_problem,problem_setting,SMALSE_params);
    
%     
%     figure(11)
%     subplot(1,2,1)
%     hold on
%     plot(cell2mat(D));
%     subplot(1,2,2)
%     hold on
%     plot(Q)
end
end
tt=tic;
fprintf('Time of calculation = %.2f seconds.\n',double(tt-t)/1e6);