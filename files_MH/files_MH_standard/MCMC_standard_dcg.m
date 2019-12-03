function [ SAMPLES,OBSERVATIONS,evaluated_u_mh,evaluated_Gu_mh, TIMES, W ] = MCMC_standard_dcg( u,y,sp,sigma,gamma,mu0,g,N,period_saving,batch,label )
%MCMC stadard algorithm without approximation
solver=@(A,b)PDCG(A,b,[],sp.W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
G=@(x)observation_2mat_solver(x,sp.M,solver);

% period_saving ... pocet snimku do ulozeni a drawnow
tic
%scatter(u(1),u(2),'*r');
cov_mat=diag(g.*g);
n=length(u);
sample=0;
[Gu,~,sol,~,~,~,freeNode]=G(u);

W=GramSchmidt_one(sp.W,sol(freeNode));
solver=@(A,b)PDCG(A,b,[],W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
G=@(x)observation_2mat_solver(x,sp.M,solver);

evaluated_u_mh=u';
evaluated_Gu_mh=Gu';
Au=(y-Gu)'*diag(1./(sigma.^2))*(y-Gu)/2;
Bu=(u-mu0)'*diag(1./(gamma.^2))*(u-mu0)/2;
prijato=0;
neprijato=0;
% rng(batch);
counter_save=0;
counter_accepted=0;
RATES=zeros(1,N/period_saving);
TIMES=zeros(1,N/period_saving);
VARS=zeros(n,N/period_saving);
SAMPLES=zeros(N,n);
OBSERVATIONS=zeros(N,length(Gu));
ITERATIONS=TIMES;
BASISSIZE=TIMES;
rng(22);
RES0=[]; ITER=[]; IMP=[];
for l=1:N
    counter_save=counter_save+1;
    p=u+chol(cov_mat)*my_normrnd(zeros(n,1),ones(n,1));
    [Gp,~,sol,~,~,~,freeNode,iter,resvec]=G(p);
%     disp(['1: res0=' num2str(resvec(1))]);
%     disp(['2: iter=' num2str(iter)]);
    if iter>1 && size(W,2)<50
        [W,imp]=GramSchmidt_one(W,sol(freeNode));
%         disp(['3: imp=' num2str(imp)]);
    end
%     RES0=[RES0; resvec(1)];
%     ITER=[ITER; iter];
%     IMP =[IMP;  imp];
    solver=@(A,b)PDCG(A,b,[],W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
    G=@(x)observation_2mat_solver(x,sp.M,solver);
    
    [lnratio,Ap,Bp,Gp]=MCMC_alpha_dcg(p,y,Gp,sigma,gamma,mu0,Au,Bu);
    evaluated_u_mh=[evaluated_u_mh; p'];      % remember all samples, for which 
    evaluated_Gu_mh=[evaluated_Gu_mh; Gp'];     % G was calculated
    if log(rand)<lnratio % p prijato
        u=p; % prijmu navrh p
        Au=Ap;
        Bu=Bp;
        Gu=Gp;
        counter_accepted=counter_accepted+1;
        if n==2
            scatter(p(1),p(2),'.k');
        else
            scatter3(p(1),p(2),p(3),'.k');
        end
        prijato=prijato+1;
        %disp(['prijato ' num2str(prijato)])
    else
        if n==2
            scatter(p(1),p(2),'.r');
        else
            scatter3(p(1),p(2),p(3),'.r');
        end
        neprijato=neprijato+1;
        %disp(['!!! neprijato ' num2str(neprijato)])
    end

    sample=sample+1;
    SAMPLES(sample,:)=u;
    OBSERVATIONS(sample,:)=Gu;
    TIMES(:,l)=toc;
    ITERATIONS(l)=iter;
    BASISSIZE(l)=size(W,2);
    if counter_save==period_saving
        counter_save=0;
        rate=counter_accepted/period_saving;
        counter_accepted=0;
%         drawnow
        RATES(:,l)=rate;
        VARS(:,l)=g;
%         [status, msg, msgID] = mkdir('res');
%         save(['res/RATES_' label '.mat'],'RATES')
%         save(['res/TIMES0_' label '.mat'],'TIMES')
%         save(['res/VARS_' label '.mat'],'VARS')
%         save(['res/SAMPLES_' label '.mat'],'SAMPLES');
%         save(['res/OBSERVATIONS_' label '.mat'],'OBSERVATIONS');
        disp([num2str(prijato) ' prijato, ' num2str(neprijato) ' neprijato']);
    end
end
figure(20); plot(TIMES); hold on
figure(21); plot(cumsum(ITERATIONS)); hold on
figure(22); plot(BASISSIZE); hold on
%plot_chain(SAMPLES,OBSERVATIONS);
%plot_chain(evaluated_u_mh,evaluated_Gu_mh);
figure(23); plot(RES0); set(gca,'YScale','log')
figure(24); plot(ITER); set(gca,'YScale','log')
figure(25); plot(IMP); set(gca,'YScale','log')
end

