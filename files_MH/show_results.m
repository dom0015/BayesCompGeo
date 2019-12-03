function [ LENG, PSTD, ARAT, NOOG, CORR, EFFI ] = show_results( sigmaMH, prijato, neprijato, nepredprijato, SAMPLES, RATES )
%SHOW_RESULTS - TABLE WITH RESULTS
% LENG = chain length
% PSTD = proposal std
% ARAT = acc rate
% NOOG = no of G evals
% CORR = corr length
% EFFI = (no of G evals) * (corr length) / (chain length)
[ending,no_par] = size(SAMPLES);
LENG = ending;
PSTD = sigmaMH(1);
ARAT = mean(RATES);
NOOG = prijato + neprijato;
% autocorrelation estimation:
starting = max(100,floor(ending/100));
CHAIN = SAMPLES(starting+1:ending,:);
temp = 0;
for j=1:no_par
    tau = autocorr_compare(CHAIN(:,j));
    temp = max(temp, tau);
end
CORR = max(temp);
EFFI = NOOG .* CORR ./ LENG;
disp([LENG' PSTD' ARAT' NOOG' CORR' EFFI'])

%figure; imagesc( [-10 -1], [-10 -1], hist2d(CHAIN(:,1:2),-10:0.1:-1,-10:0.1:-1));
%figure; plot(CHAIN(:,1),CHAIN(:,2),'.k');
end

