function [Q,D,PRESSURE,ugrad] = coupled_solver(params,hydro_problem,problem_setting,SMALSE_params)


SMALSE_params.rel=1.0e-2;
hydro_problem.mat_frac=params;
%% depends on fracture aperture
%d = exp([-6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7 -6 -8.5 -7])/100;
no_fractures=hydro_problem.no_fractures;
lengths=hydro_problem.lengths;
d = 1e-4*ones(no_fractures,1);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end
D_old=D;
eps_coupling=SMALSE_params.eps_coupling;
res_d=[];
res_press=[];
alpha=0.5;
beta=1e-1;
tic;
rel_D_=[];
rel_P_=[];
for i=1:SMALSE_params.coupling_iter    
    [PRESSURE,~,ugrad,Q,PRESSURE_diff]=tocouple_handle(D,hydro_problem);
    res_press(:,i)=cell2mat(PRESSURE);
    
    if i>1
        for j=1:size(PRESSURE_diff,1)
            for jj=1:size(PRESSURE_diff,2)
                PRESSURE_diff{j,jj}=(PRESSURE_diff{j,jj}*alpha+(1-alpha)*PRESSURE_old{j,jj});
            end 
        end
        ugrad=ugrad*alpha+(1-alpha)*ugrad_old;
    end
    for j=1:size(PRESSURE_diff,1)
        for jj=1:size(PRESSURE_diff,2)
            PRESSURE_diff{j,jj}=(PRESSURE_diff{j,jj}*beta);
        end 
    end
    ugrad=beta*ugrad;
    [problem_setting] = assembly_FETI_frac_rhs(problem_setting,PRESSURE_diff,-1*ugrad);
    
 
    [D,problem_setting] = SMALSE_solver(problem_setting,SMALSE_params);
    
    
    for j=1:length(D)
        D{j}=max(D{j},1e-10);
    end
    
    res_d(:,i)=cell2mat(D);
    if i>1
        for j=1:length(D)
            D{j}=(D{j}*alpha+(1-alpha)*D_old{j});
        end
    end
    
    tmp=sqrt(sum((res_press(:,1:end-1)-res_press(:,2:end)).^2,1))./sqrt(sum((res_press(:,2:end)).^2,1));
    tmp2=sqrt(sum((res_d(:,1:end-1)-res_d(:,2:end)).^2,1))./sqrt(sum((res_d(:,2:end)).^2,1));
    
    D_old=D;
    PRESSURE_old=PRESSURE_diff;
    ugrad_old=ugrad;
    
    
    % SMALSE precision to match coupled precision
    if i>2
        SMALSE_params.rel=tmp2(end)/1000;
    end
    % coupling parameters update
    alpha=min(0.8,alpha^(0.5));
%     disp(["alpha", num2str(alpha)]);
    beta=min(1,beta^(0.5)+1e-1);
    if i>1
        if beta==1 && alpha>=0.1 && tmp(end)<eps_coupling && tmp2(end)<eps_coupling
            break
        end
    end
    
    if beta==1
        rel_D_=[rel_D_ tmp(end)];
        rel_P_=[rel_P_ tmp2(end)];
    end
end

%     figure(5)
%     subplot(1,2,1)
%     hold on
%     plot(rel_P_,'k-','LineWidth',1)
%     set(gca,'YScale','log')
%     %figure(5)
%     subplot(1,2,2)
%     yyaxis right
%     hold on
%     plot(rel_D_,'k-','LineWidth',1)
%     set(gca,'YScale','log')


    if i>1
%         fprintf('i=%d: alpha = %.2f , beta = %.2f || pres = %d || d = %d\n',i,alpha,beta,tmp(end),tmp2(end));
    end
end

