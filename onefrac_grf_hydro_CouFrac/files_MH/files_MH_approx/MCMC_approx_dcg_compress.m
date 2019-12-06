function [ SAMPLES,OBSERVATIONS,MULTIPLICITY, G_approx, log_approx, evaluated_Gu_mh, TIMES, W, CUM_ITER, BASISSIZE ] = MCMC_approx_dcg_compress( u,y,sp,sigma,gamma,mu0,g,N,period_saving,batch,label,evaluated_u_init,evaluated_Gu_init,evaluated_M_init,upgrade_approx,approx_type )
%MCMC using an approximation
%   uses an approximation of the observation operator G
%   to reduce the time consumption
% scatter(u(1),u(2),'*r'); hold on
if isstruct(sp)
    flag_iterative=1;
else
    flag_iterative=0;
end
if flag_iterative
    threshold_iter=4;
    threshold_size=100;
    solver=@(A,b)PDCG(A,b,[],sp.W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
    G=@(x)observation_2mat_solver(x,sp.M,solver);
else
    G=sp;
    W=[];
end

tic
cov_mat=diag(g.*g);
n=length(u);

if flag_iterative
    [Gu,~,sol,~,~,~,freeNode,iter,resvec]=G(u);
    W=sp.W;
    if iter>threshold_iter && size(W,2)<threshold_size
        [W,imp]=GramSchmidt_one(W,sol(freeNode));
    end

    solver=@(A,b)PDCG(A,b,[],W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
    G=@(x)observation_2mat_solver(x,sp.M,solver);
else
    Gu=G(u);
end

evaluated_u_mh=u';
evaluated_Gu_mh=Gu';
%% initial approximation based on initial MH run
if approx_type==0
    G_approx=upgrade_approx(evaluated_u_init,evaluated_Gu_init); % RBF
elseif approx_type==1
    G_approx=upgrade_approx(evaluated_u_init,evaluated_Gu_init,evaluated_M_init); % for COL
else
    WW=W;
    if size(W,2)>20
        WW=W(:,end-20+1:end);
    end
    solver_=@(A,b)PDCG(A,b,[],WW,[],sp.prec,sp.cg_accuracy,0);
    upgrade_approx=@(a,b)(@(x)observation_2mat_solver(x,sp.M,solver_));
    G_approx=upgrade_approx();
end
G_approx_u=G_approx(u);
Au_approx=(y-G_approx_u)'*diag(1./(sigma.^2))*(y-G_approx_u)/2;
Au=(y-Gu)'*diag(1./(sigma.^2))*(y-Gu)/2;
Bu=(u-mu0)'*diag(1./(gamma.^2))*(u-mu0)/2;
prijato=0;
neprijato=0;
nepredprijato=0;
counter_save=0;
counter_accepted=0;
log_approx=[];
RATES=zeros(1,N/period_saving);
TIMES=zeros(1,N);
ITERATIONS=TIMES;
BASISSIZE=TIMES;
VARS=zeros(n,N/period_saving);
recent_sample=1;
SAMPLES=zeros(period_saving,n); SAMPLES(recent_sample,:)=u';
MULTIPLICITY=zeros(period_saving,1); %MULTIPLICITY(recent_sample)=1;
OBSERVATIONS=zeros(period_saving,length(Gu)); OBSERVATIONS(recent_sample,:)=Gu';
rng(batch)
for l=1:N
    iter=0;
    counter_save=counter_save+1;
    p=u+chol(cov_mat)*my_normrnd(zeros(n,1),ones(n,1));
    [lnratio_approx,Ap_approx,Bp,A_approx,Gp_app]=MCMC_alpha_approx(p,y,G_approx,sigma,gamma,mu0,Au_approx,Bu);
    if log(rand)<lnratio_approx % predprijato
        if flag_iterative
            [Gp,~,sol,~,~,~,freeNode,iter]=G(p);
            if iter>threshold_iter && size(W,2)<threshold_size
    %             disp(iter)
                W=GramSchmidt_one(W,sol(freeNode));
            end
            solver=@(A,b)PDCG(A,b,[],W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
            G=@(x)observation_2mat_solver(x,sp.M,solver);
        else
            Gp=G(p);
        end
        
        [lnratio,Ap,Gp]=MCMC_alpha_exact_dcg(y,Gp,sigma,Au,A_approx);
        %disp([Gp_app'; Gp'])
        if approx_type==0 % RBF
            evaluated_u_mh=[evaluated_u_mh; p'];      % remember all samples, for which 
            evaluated_Gu_mh=[evaluated_Gu_mh; Gp'];     % G was calculated
            ff=1;
            evaluated_u=[evaluated_u_init; evaluated_u_mh(1:ff:end,:)];
            evaluated_Gu=[evaluated_Gu_init; evaluated_Gu_mh(1:ff:end,:)];
            G_approx=upgrade_approx(evaluated_u,evaluated_Gu);
        elseif approx_type==1 % COL
            evaluated_u=[evaluated_u_init; SAMPLES(1:recent_sample,:); p'];
            evaluated_Gu=[evaluated_Gu_init; OBSERVATIONS(1:recent_sample,:); Gp'];
            evaluated_M=[evaluated_M_init; MULTIPLICITY(1:recent_sample,:); 1];
            G_approx=upgrade_approx(evaluated_u,evaluated_Gu,evaluated_M);
        else
            WW=W;
            if size(W,2)>20
                WW=W(:,end-20+1:end);
            end
            solver_=@(A,b)PDCG(A,b,[],WW,[],sp.prec,sp.cg_accuracy,0);
            upgrade_approx=@(a,b)(@(x)observation_2mat_solver(x,sp.M,solver_));
            G_approx=upgrade_approx();
        end
        if log(rand)<lnratio % p prijato
            u=p; % prijmu navrh p
            Au_approx=Ap_approx;
            Au=Ap;
            Bu=Bp;
            Gu=Gp;
            recent_sample=recent_sample+1;
            MULTIPLICITY(recent_sample)=1;
            SAMPLES(recent_sample,:)=u;
            OBSERVATIONS(recent_sample,:)=Gu;
            counter_accepted=counter_accepted+1;
%             if n==2
%                 scatter(p(1),p(2),'.k');
%             else
%                 scatter3(p(1),p(2),p(3),'.k');
%             end
            prijato=prijato+1;
            %disp(['prijato ' num2str(prijato)])
        else
            MULTIPLICITY(recent_sample)=MULTIPLICITY(recent_sample)+1;
%             if n==2
%                 scatter(p(1),p(2),'.r');
%             else
%                 scatter3(p(1),p(2),p(3),'.r');
%             end
            neprijato=neprijato+1;
            %disp([Gp Gp_app]);
            %disp(['!!! neprijato ' num2str(neprijato)])
            log_approx=[log_approx; l];
        end
    else
        MULTIPLICITY(recent_sample)=MULTIPLICITY(recent_sample)+1;
%         if n==2
%             scatter(p(1),p(2),'.g');
%         else
%             scatter3(p(1),p(2),p(3),'.g');
%         end
        nepredprijato=nepredprijato+1;
        %disp(['nepredprijato ' num2str(nepredprijato) ' ' num2str(lnratio)])
    end
    TIMES(l)=toc;
    ITERATIONS(l)=iter;
    BASISSIZE(l)=size(W,2);
    
    if counter_save==period_saving
        counter_save=0;
        rate=counter_accepted/period_saving;
        counter_accepted=0;
        %drawnow
        RATES(:,l/period_saving)=rate;
        VARS(:,l/period_saving)=g;
        SAMPLES=SAMPLES(1:recent_sample,:);
        OBSERVATIONS=OBSERVATIONS(1:recent_sample,:);
        MULTIPLICITY=MULTIPLICITY(1:recent_sample,:);
        [status, msg, msgID] = mkdir('res');
        save(['res/DAMH_' label '.mat'],'SAMPLES','OBSERVATIONS',...
            'MULTIPLICITY','RATES','VARS','TIMES',...
            'prijato','neprijato','nepredprijato','-v7.3');
        SAMPLES=[SAMPLES; zeros(period_saving,n)];
        MULTIPLICITY=[MULTIPLICITY; zeros(period_saving,1)];
        OBSERVATIONS=[OBSERVATIONS; zeros(period_saving,length(Gu))];
        disp([num2str(prijato) ' prijato, ' num2str(neprijato) ' neprijato, ' num2str(nepredprijato) ' nepredprijato']);
    end
end
%plot_chain(SAMPLES,OBSERVATIONS);
[ LENG, PSTD, ARAT, NOOG, CORR, EFFI ] = show_results( g, prijato, neprijato, nepredprijato, SAMPLES, RATES );
save(['res/INFO_' label '.mat'],'g','neprijato','prijato','nepredprijato','LENG','PSTD','ARAT','NOOG','CORR','EFFI')
% figure; imagesc(W)
%a=min(SAMPLES(:,1)); b=max(SAMPLES(:,1));
%c=min(SAMPLES(:,1)); d=max(SAMPLES(:,1));
%figure; imagesc(hist2d(SAMPLES,linspace(c,d,100),linspace(a,b,100)))

figure(23); plot(TIMES); title('times'); hold on
CUM_ITER=cumsum(ITERATIONS);
figure(24); plot(CUM_ITER); title('iter'); hold on
figure(25); plot(BASISSIZE); title('basis size'); hold on
end        

