function [Q_neg,additional_output] = wrapper_hydro(log_mat_frac,param,hydro_problem,data_generator)
%WRAPPER_HYDRO Summary of this function goes here
%   Detailed explanation goes here
mat_frac=10^log_mat_frac;
hydro_problem.mat_frac=mat_frac;



D{1}=10.^(data_generator.gen_basis*param+data_generator.prior_mean);

hold on
plot(D{1})
drawnow

[PRESSURE,u0,GRAD,Q,PRESSURE_diff] = tocouple_handle(D,hydro_problem);
if length(Q)>4
    Q_neg=Q(2:end)';
else
    Q_neg=Q';
end
additional_output=[];
end

