function [] = viz_corrplot_all(batches,label,u_real,gamma,mu0,folder_name)
%VIZ_CORRPLOT Summary of this function goes here
%   Detailed explanation goes here


MH_SAMPLES=[];
MH_MULTIPLICITY=[];
MH_u=[];
MH_Gu=[];
DAMHSMU_SAMPLES=[];
DAMHSMU_u=[];
DAMHSMU_Gu=[];
DAMH_SAMPLES=[];
DAMH_MULTIPLICITY=[];
DAMH_u=[];
DAMH_Gu=[];
for i=1:batches
    filename = choose0or1('MH', i, folder_name, label);
    load(filename)
    MH_SAMPLES=[MH_SAMPLES; SAMPLES];
    MH_MULTIPLICITY=[MH_MULTIPLICITY; MULTIPLICITY];
    MH_u=[MH_u; evaluated_u_mh];
    MH_Gu=[MH_Gu; evaluated_Gu_mh];
    filename = choose0or1('DAMHSMUdata', i, folder_name, label);
    load(filename)
%     DAMHSMU_SAMPLES=[DAMHSMU_SAMPLES; SAMPLES];
    DAMHSMU_u=[DAMHSMU_u; evaluated_u_mh];
    DAMHSMU_Gu=[DAMHSMU_Gu; evaluated_Gu_mh];
    filename = choose0or1('DAMHdata', i, folder_name, label);
    load(filename)
%     DAMH_SAMPLES=[DAMH_SAMPLES; SAMPLES];
    DAMH_u=[DAMH_u; evaluated_u_mh];
    DAMH_Gu=[DAMH_Gu; evaluated_Gu_mh];
    filename = choose0or1('DAMH', i, folder_name, label);
    load(filename)
    DAMH_SAMPLES=[DAMH_SAMPLES; SAMPLES];
    DAMH_MULTIPLICITY=[DAMH_MULTIPLICITY; MULTIPLICITY];
end
ALL_u=[MH_u; DAMHSMU_u; DAMH_u];
% ALL_Gu=[MH_Gu; DAMHSMU_Gu; DAMH_Gu];
prior_std=gamma;
prior_mean=mu0;
save(['allres/DAMH_ALL_' label '.mat'],'DAMH_SAMPLES','DAMH_MULTIPLICITY','u_real','prior_std','prior_mean')
save(['allres/MH_ALL_' label '.mat'],'MH_SAMPLES','MH_MULTIPLICITY','u_real','prior_std','prior_mean')

figure
corrplot(MH_u);
figure
corrplot(DAMHSMU_u);
figure
corrplot(DAMH_u);
figure
corrplot(MH_SAMPLES);
figure
corrplot(DAMH_SAMPLES);