function [ G_approx, log_approx, evaluated_Gu_mh ] = MCMC_approx( u,y,G,sigma,gamma,mu0,g,N,period_saving,batch,label,evaluated_u_init,evaluated_Gu_init,upgrade_approx,approx_type )
%MCMC using an approximation
%   uses an approximation of the observation operator G
%   to reduce the time consumption
% scatter(u(1),u(2),'*r'); hold on
tic
cov_mat=diag(g.*g);
n=length(u);
sample=0;
Gu=G(u);
evaluated_u_mh=u';
evaluated_Gu_mh=Gu';
%% initial approximation based on initial MH run
G_approx=upgrade_approx(evaluated_u_init,evaluated_Gu_init); % for COL
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
VARS=zeros(n,N/period_saving);
SAMPLES=zeros(N,n);
OBSERVATIONS=zeros(N,length(Gu));
rng(batch)
for l=1:N
    counter_save=counter_save+1;
    p=u+chol(cov_mat)*my_normrnd(zeros(n,1),ones(n,1));
    [lnratio_approx,Ap_approx,Bp,A_approx,Gp_app]=MCMC_alpha_approx(p,y,G_approx,sigma,gamma,mu0,Au_approx,Bu);
    if log(rand)<lnratio_approx % predprijato
        [lnratio,Ap,Gp]=MCMC_alpha_exact(p,y,G,sigma,Au,A_approx);
        %disp([Gp_app'; Gp'])
        evaluated_u_mh=[evaluated_u_mh; p'];      % remember all samples, for which 
        evaluated_Gu_mh=[evaluated_Gu_mh; Gp'];     % G was calculated
        if approx_type==0 % RBF
            ff=1;
            evaluated_u=[evaluated_u_init; evaluated_u_mh(1:ff:end,:)];
            evaluated_Gu=[evaluated_Gu_init; evaluated_Gu_mh(1:ff:end,:)];
            G_approx=upgrade_approx(evaluated_u,evaluated_Gu);
        else % COL
            evaluated_u=[evaluated_u_init; SAMPLES(1:sample,:); p'];
            evaluated_Gu=[evaluated_Gu_init; OBSERVATIONS(1:sample,:); Gp'];
            G_approx=upgrade_approx(evaluated_u,evaluated_Gu);
        end
        if log(rand)<lnratio % p prijato
            u=p; % prijmu navrh p
            Au_approx=Ap_approx;
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
            %disp([Gp Gp_app]);
            %disp(['!!! neprijato ' num2str(neprijato)])
            log_approx=[log_approx; l];
        end
    else
        if n==2
            scatter(p(1),p(2),'.g');
        else
            scatter3(p(1),p(2),p(3),'.g');
        end
        nepredprijato=nepredprijato+1;
        %disp(['nepredprijato ' num2str(nepredprijato) ' ' num2str(lnratio)])
    end
    sample=sample+1;
    SAMPLES(sample,:)=u;
    OBSERVATIONS(sample,:)=Gu;
    TIMES(l)=toc;
    
    if counter_save==period_saving
        counter_save=0;
        rate=counter_accepted/period_saving;
        counter_accepted=0;
        %drawnow
        RATES(:,l/period_saving)=rate;
        VARS(:,l/period_saving)=g;
%         [status, msg, msgID] = mkdir('res');
%         SAMPLES_save = SAMPLES(l-period_saving+1:l,:);
%         OBSERVATIONS_save = OBSERVATIONS(l-period_saving+1:l,:);
%         save(['res/RATES_' label '.mat'],'RATES','-v7.3')
%         save(['res/TIMES_' label '.mat'],'TIMES','-v7.3')
%         save(['res/VARS_' label '.mat'],'VARS','-v7.3')
%         save(['res/SAMPLES_' num2str(l/period_saving) '_' label '.mat'],'SAMPLES_save','-v7.3');
%         save(['res/OBSERVATIONS_' num2str(l/period_saving) '_' label '.mat'],'OBSERVATIONS_save','-v7.3');
        disp([num2str(prijato) ' prijato, ' num2str(neprijato) ' neprijato, ' num2str(nepredprijato) ' nepredprijato']);
    end
end

figure(23); plot(TIMES); hold on
%plot_chain(SAMPLES,OBSERVATIONS);
[ LENG, PSTD, ARAT, NOOG, CORR, EFFI ] = show_results( g, prijato, neprijato, nepredprijato, SAMPLES, RATES );
save(['res/INFO_' label '.mat'],'g','neprijato','prijato','nepredprijato','LENG','PSTD','ARAT','NOOG','CORR','EFFI')

end        

