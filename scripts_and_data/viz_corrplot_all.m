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
    cd([folder_name '/res'])
    t0=dir(['MH_' label '_' num2str(i) '_0.mat']);
    t1=dir(['MH_' label '_' num2str(i) '_1.mat']);
    cd('..');
    cd('..');
    if t0.bytes > t1.bytes 
        filename=[folder_name '/res/MH_' label '_' num2str(i) '_0.mat'];
        disp('MH0')
    else
        filename=[folder_name '/res/MH_' label '_' num2str(i) '_1.mat'];
        disp('MH1')
    end
    
    load(filename)
    MH_SAMPLES=[MH_SAMPLES; SAMPLES];
    MH_MULTIPLICITY=[MH_MULTIPLICITY; MULTIPLICITY];
    MH_u=[MH_u; evaluated_u_mh];
    MH_Gu=[MH_Gu; evaluated_Gu_mh];
    if t0.bytes > t1.bytes 
        filename=[folder_name '/res/DAMHSMUdata_' label '_' num2str(i) '_0.mat'];
        disp('MH0')
    else
        filename=[folder_name '/res/DAMHSMUdata_' label '_' num2str(i) '_1.mat'];
        disp('MH1')
    end
    load(filename)
%     DAMHSMU_SAMPLES=[DAMHSMU_SAMPLES; SAMPLES];
    DAMHSMU_u=[DAMHSMU_u; evaluated_u_mh];
    DAMHSMU_Gu=[DAMHSMU_Gu; evaluated_Gu_mh];
    if t0.bytes > t1.bytes 
        filename=[folder_name '/res/DAMHdata_' label '_' num2str(i) '_0.mat'];
        disp('MH0')
    else
        filename=[folder_name '/res/DAMHdata_' label '_' num2str(i) '_1.mat'];
        disp('MH1')
    end
    load(filename)
%     DAMH_SAMPLES=[DAMH_SAMPLES; SAMPLES];
    DAMH_u=[DAMH_u; evaluated_u_mh];
    DAMH_Gu=[DAMH_Gu; evaluated_Gu_mh];
    cd([folder_name '/res'])
    t0=dir(['DAMH_' label '_' num2str(i) '_0.mat']);
    t1=dir(['DAMH_' label '_' num2str(i) '_1.mat']);
    cd('..');
    cd('..');
    if t0.bytes>t1.bytes
        filename=[folder_name '/res/DAMH_' label '_' num2str(i) '_0.mat'];
        disp('DAMMH0')
    else
        filename=[folder_name '/res/DAMH_' label '_' num2str(i) '_1.mat'];
        disp('DAMMH1')
    end
    
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