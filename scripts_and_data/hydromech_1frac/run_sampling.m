% parallel.importProfile('/apps/all/MATLAB/2015b-EDU/SalomonPBSPro.settings')
% cluster = parcluster('local');
% pool=parpool(cluster,16);
%parpool(4)
%% Bayesian solution of an inverse problem
addpath(genpath('../../files_MH'))
addpath(genpath('../../files_elasticity'))
addpath(genpath('../../files_hydro'))
addpath(genpath('../../files_other'))

Nbatch=15;

ulohy=cell(1,1);
ulohy{1}=load('uloha2.mat');
ulohy{1}.SMALSE_params.coupling_iter=100;

ulohy{1}.SMALSE_params.eps_coupling=1e-6;
ulohy{1}.SMALSE_params.print_couple=true;
ulohy{1}.SMALSE_params.print=true;
G=@(x)G_operator(10.^x,ulohy);

ALL_D_cell=cell(Nbatch,1);
ALL_Q_cell=cell(Nbatch,1);
param=linspace(-4.5,-5,31);
param=param(1:30);
param=reshape(param,2,15);

for batch=1:Nbatch
    Q_temp=[];
    D_temp=[];
    param_temp=param(:,batch)';
    for i=1:2
        param_temp(i)
        [Q,D]=G(param_temp(i));
        Q_temp(i,:)=Q';
        D_temp(i,:)=D';
    end
    ALL_Q_cell{batch}=Q_temp;
    ALL_D_cell{batch}=D_temp;
end

ALL_Q=cell2mat(ALL_Q_cell);
ALL_D=cell2mat(ALL_D_cell);

%save('hydromech_1frac.mat','ALL_Q','ALL_D')