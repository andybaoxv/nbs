% This script apply NBS

% Set the number of latent factors
option.K = 4;
option.iter = 500;

% Defaule clustering settings
option.zoptions.iter = 1000;

% Relate the clustering performance
option.zoptions.gamma = 200;
option.zoptions.tof = 1e-4;
option.zoptions.dis = false;
option.zoptions.distance = 'nnls';

tic
[W,H,tstats] =  nbs_nmf_cluster_network_nmf(sub_gene_data,option.K,0,...
    -1,option.zoptions.gamma,sparse(subnet_laplacian),option.zoptions);
toc
indClust = NMFCluster(H);

% Compute cluster size
counter_NBS = zeros(1,4);
for i = 1:length(indClust)
    counter_NBS(indClust(i)) = counter_NBS(indClust(i))+1; 
end

counter_NBS
tstats