function [ SAMPLES,OBSERVATIONS,MULTIPLICITY,evaluated_u_mh,evaluated_Gu_mh, TIMES, W, CUM_ITER, BASISSIZE ] = MCMC_standard_dcg_compress( u,y,sp,sigma,gamma,mu0,g,N,period_saving,batch,label )
%MCMC stadard algorithm without approximation
if isstruct(sp)
    flag_iterative=1;
else
    flag_iterative=0;
end
if flag_iterative
    threshold_iter=6;
    threshold_size=50;
    solver=@(A,b)PDCG(A,b,[],sp.W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
    G=@(x)observation_2mat_solver(x,sp.M,solver);
else
    G=sp;
    W=[];
end

% period_saving ... pocet snimku do ulozeni a drawnow
tic
%scatter(u(1),u(2),'*r');
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
Au=(y-Gu)'*diag(1./(sigma.^2))*(y-Gu)/2;
Bu=(u-mu0)'*diag(1./(gamma.^2))*(u-mu0)/2;
prijato=0;
neprijato=0;
% rng(batch);
counter_save=0;
counter_accepted=0;
RATES=zeros(1,N/period_saving);
TIMES=zeros(1,N);
VARS=zeros(n,N/period_saving);
recent_sample=1;
SAMPLES=zeros(period_saving,n); SAMPLES(recent_sample,:)=u';
MULTIPLICITY=zeros(period_saving,1); MULTIPLICITY(recent_sample)=1;
OBSERVATIONS=zeros(period_saving,length(Gu)); OBSERVATIONS(recent_sample,:)=Gu';
ITERATIONS=TIMES;
BASISSIZE=TIMES;
rng(22);
for l=1:N
    counter_save=counter_save+1;
    p=u+chol(cov_mat)*my_normrnd(zeros(n,1),ones(n,1));
    if flag_iterative
        [Gp,~,sol,~,~,~,freeNode,iter,resvec]=G(p);
    %     disp(['1: res0=' num2str(resvec(1))]);
    %     disp(['2: iter=' num2str(iter)]);
        if iter>threshold_iter && size(W,2)<threshold_size
            [W,imp]=GramSchmidt_one(W,sol(freeNode));
    %         disp(['3: imp=' num2str(imp)]);
        end
    %     RES0=[RES0; resvec(1)];
    %     ITER=[ITER; iter];
    %     IMP =[IMP;  imp];
        solver=@(A,b)PDCG(A,b,[],W,[],sp.prec,sp.cg_accuracy,sp.max_iter);
        G=@(x)observation_2mat_solver(x,sp.M,solver);
    else
        Gp=G(p);
    end
    
    [lnratio,Ap,Bp,Gp]=MCMC_alpha_dcg(p,y,Gp,sigma,gamma,mu0,Au,Bu);
%     evaluated_u_mh=[evaluated_u_mh; p'];      % remember all samples, for which 
%     evaluated_Gu_mh=[evaluated_Gu_mh; Gp'];     % G was calculated
    if log(rand)<lnratio % p prijato
        u=p; % prijmu navrh p
        Au=Ap;
        Bu=Bp;
        Gu=Gp;
        recent_sample=recent_sample+1;
        MULTIPLICITY(recent_sample)=1;
        SAMPLES(recent_sample,:)=u;
        OBSERVATIONS(recent_sample,:)=Gu;
        counter_accepted=counter_accepted+1;
        prijato=prijato+1; 
%         if n==2
%             scatter(p(1),p(2),'.k');
%         else
%             scatter3(p(1),p(2),p(3),'.k');
%         end
        %disp(['prijato ' num2str(prijato)])
    else
%         if n==2
%             scatter(p(1),p(2),'.r');
%         else
%             scatter3(p(1),p(2),p(3),'.r');
%         end
        MULTIPLICITY(recent_sample)=MULTIPLICITY(recent_sample)+1;
        neprijato=neprijato+1;
        %disp(['!!! neprijato ' num2str(neprijato)])
    end

    TIMES(:,l)=toc;
    if flag_iterative
        ITERATIONS(l)=iter;
        BASISSIZE(l)=size(W,2); disp([size(W,2)  iter])
    end
    if counter_save==period_saving
        counter_save=0;
        rate=counter_accepted/period_saving;
        counter_accepted=0;
%         drawnow
        RATES(:,l/period_saving)=rate;
        VARS(:,l/period_saving)=g;
        SAMPLES=SAMPLES(1:recent_sample,:);
        OBSERVATIONS=OBSERVATIONS(1:recent_sample,:);
        MULTIPLICITY=MULTIPLICITY(1:recent_sample,:);
        [status, msg, msgID] = mkdir('res');
        save(['res/MH_' label '.mat'],'SAMPLES','OBSERVATIONS',...
            'MULTIPLICITY','RATES','VARS','TIMES',...
            'prijato','neprijato','-v7.3');
        SAMPLES=[SAMPLES; zeros(period_saving,n)];
        MULTIPLICITY=[MULTIPLICITY; zeros(period_saving,1)];
        OBSERVATIONS=[OBSERVATIONS; zeros(period_saving,length(Gu))];
        disp([num2str(prijato) ' prijato, ' num2str(neprijato) ' neprijato']);
    end
end

figure(20); plot(TIMES); title('times'); hold on
if flag_iterative
    CUM_ITER=cumsum(ITERATIONS);
    figure(21); plot(CUM_ITER); title('iter'); hold on
    figure(22); plot(BASISSIZE); title('basis size'); hold on
end
%plot_chain(SAMPLES,OBSERVATIONS);
%plot_chain(evaluated_u_mh,evaluated_Gu_mh);
end

