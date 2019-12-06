%load('chain_example.mat'); X=SAMPLES(:,1);

N=1e4;
X=zeros(N,1);
for i=2:N
    prop=X(i-1)+normrnd(0,1);
    if rand>(abs(prop)/10)
        X(i)=prop;
    else
        X(i)=X(i-1);
    end
end
figure; plot(X)
L=length(X);

%% Batch means - A Comparison of Methods for Computing Autocorrelation Time, Madeleine B. Thompson
s2=var(X);
batch_size=floor(L^(2/3));
no_batches=floor(L^(1/3));
temp=[];
for j=1:batch_size:(L-batch_size)
    batch_mean=mean(X(j:(j+batch_size-1)));
    temp=[temp batch_mean];
end
s2m=var(temp);
tau_batchmeans=batch_size*s2m/s2;

%% Autocorrelation time estimation, Dan Foreman-Mackey
M=floor(L/100);
cf_matlab=autocorr(X,M);
tau_Foreman=1+2*sum(cf_matlab);

%% Stable estimates of autocorrelation length from short GSP-Phot Aeneas MCMC chains
XX=X-mean(X);
temp1=XX(2:end).*XX(1:end-1);
temp2=XX(1:end-1).^2;
phi=sum(temp1)/sum(temp2);
tau_stable_exp=-2*1/log(phi);
tau_stable_int=2*1/(1-phi)-2*1/2;

% %% The Yule Walker Equations for the AR Coefficients, Gidon Eshel
% M=min(floor(L/10),10000);
% cf_matlab=autocorr(X,M);
% R=toeplitz(cf_matlab(1:M));
% rhs=cf_matlab(2:end);
% pi_AR=R\rhs;
% tau_AR=(1-rhs'*pi_AR)/( (1-ones(1,M)*pi_AR).^2 );

tau=[tau_batchmeans, tau_Foreman, tau_stable_exp, tau_stable_int];%, tau_AR];
%tau=tau_stable_exp;
disp(tau)

n=900;
a=zeros(1,n);
for i=1:n
    X1=X(1:end-i);
    X2=X(1+i:end);
    a(i)=corr(X1,X2);
end
figure; plot(a)