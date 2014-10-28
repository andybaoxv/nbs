% This script plots dispersion coefficient

range_factors = [2,3,4,5];

dc_NBS = [0.9761,0.9780,0.9490,0.7510];
dc_NMF = [0.9627,0.6910,0.5775,0.4609];

figure(1)
plot(range_factors,dc_NBS,'-ro',range_factors,dc_NMF,'-.b');
legend('NBS Score','NMF Score');
xlabel('Number of Latent Factors');
ylabel('Stability of Clustering Solution');
