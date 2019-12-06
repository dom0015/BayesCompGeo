%% TABLE WITH RESULTS
% LENG = chain length
% PSTD = proposal std
% ARAT = acc rate
% NOOG = no of G evals
% CORR = corr length
% EFFI = (no of G evals) * (corr length) / (chain length)
clear all
label = 'aom10e8';
no_par = 7;

for i=1:24
    i
    load(['RATES_' label '_' num2str(i) '_s.mat' ])
    load(['SAMPLES_' label '_' num2str(i) '_s.mat' ])
    load(['STATUS_' label '_' num2str(i) '_s.mat' ])
    ending = l*M;
    LENG(i) = ending;
    PSTD(i) = g(1);
    ARAT(i) = mean(RATES);
    NOOG(i) = prijato + neprijato;
    % autocorrelation estimation:
    starting = max(100,floor(ending/100));
    CHAIN = SAMPLES(starting+1:ending,:);
    temp = 0;
    for j=1:no_par
        tau = autocorr_compare(CHAIN(:,j));
        temp = max(temp, tau)
    end
    CORR(i) = max(temp);
%     plot(SAMPLES)
%     drawnow
%     pause
end

EFFI = NOOG .* CORR ./ LENG;

disp([LENG' PSTD' ARAT' NOOG' CORR' EFFI'])

plot(PSTD,EFFI,'*')
    