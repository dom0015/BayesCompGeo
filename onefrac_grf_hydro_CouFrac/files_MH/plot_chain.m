function [ ] = plot_chain( SAMPLES, OBSERVATIONS)
S=cumsum(SAMPLES);
N=(1:length(S))';
figure; plot(OBSERVATIONS);
title('observations')
figure; plot(SAMPLES)
title('samples and averages')
hold on
plot(S./repmat(N,1,size(S,2)))