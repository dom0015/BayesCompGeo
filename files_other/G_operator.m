function [observations,aperatures] = G_operator(params,ulohy)

no_confs=length(ulohy);
observations=cell(no_confs,1);
aperatures=cell(no_confs,1);
for c_ulohy=1:no_confs
    problem_setting=ulohy{c_ulohy}.problem_setting;
    hydro_problem=ulohy{c_ulohy}.hydro_problem;
    SMALSE_params=ulohy{c_ulohy}.SMALSE_params;
    SMALSE_params.eps_coupling=1e-3;    
    
    [Q,D] = coupled_solver(params,hydro_problem,problem_setting,SMALSE_params);
   observations{c_ulohy}=Q(2:end)';
   aperatures{c_ulohy}=cell2mat(D);
end
observations=cell2mat(observations);
aperatures=cell2mat(aperatures);
end

