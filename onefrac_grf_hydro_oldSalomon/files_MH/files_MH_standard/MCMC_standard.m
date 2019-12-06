function [ SAMPLES,OBSERVATIONS,evaluated_u_mh,evaluated_Gu_mh, TIMES ] = MCMC_standard( u,y,G,sigma,gamma,mu0,g,N,period_saving,batch,label )
%MCMC stadard algorithm without approximation
% period_saving ... pocet snimku do ulozeni a drawnow
tic
%scatter(u(1),u(2),'*r');
cov_mat=diag(g.*g);
n=length(u);
sample=0;
Gu=G(u);
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
rng(22);
for l=1:N
    counter_save=counter_save+1;
    p=u+chol(cov_mat)*my_normrnd(zeros(n,1),ones(n,1));
    [lnratio,Ap,Bp,Gp]=MCMC_alpha(p,y,G,sigma,gamma,mu0,Au,Bu);
    evaluated_u_mh=[evaluated_u_mh; p'];      % remember all samples, for which 
    evaluated_Gu_mh=[evaluated_Gu_mh; Gp'];     % G was calculated
    if log(rand)<lnratio % p prijato
        u=p; % prijmu navrh p
        Au=Ap;
        Bu=Bp;
        Gu=Gp;
        counter_accepted=counter_accepted+1;
%        scatter(p(1),p(2),'.k');
        prijato=prijato+1;
        %disp(['prijato ' num2str(prijato)])
    else
%        scatter(p(1),p(2),'.r');
        neprijato=neprijato+1;
        %disp(['!!! neprijato ' num2str(neprijato)])
    end

    sample=sample+1;
    SAMPLES(sample,:)=u;
    OBSERVATIONS(sample,:)=Gu;
    TIMES(:,l)=toc;
    if counter_save==period_saving
        counter_save=0;
        rate=counter_accepted/period_saving;
        counter_accepted=0;
%         drawnow
        RATES(:,l)=rate;
        VARS(:,l)=g;
        [status, msg, msgID] = mkdir('res');
        save(['res/RATES_' label '.mat'],'RATES')
        save(['res/TIMES0_' label '.mat'],'TIMES')
        save(['res/VARS_' label '.mat'],'VARS')
        save(['res/SAMPLES_' label '.mat'],'SAMPLES');
        save(['res/OBSERVATIONS_' label '.mat'],'OBSERVATIONS');
        disp([num2str(prijato) ' prijato, ' num2str(neprijato) ' neprijato']);
    end
end

%plot_chain(SAMPLES,OBSERVATIONS);
%plot_chain(evaluated_u_mh,evaluated_Gu_mh);

end

