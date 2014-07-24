% This script run consensus clustering with NBS

n_consensus = 100;
flag_normalization = true;

%% Parameter setting for NBS
% Set the number of latent factors
option.K = 4;
option.iter = 500;

% Defaule clustering settings
option.zoptions.iter = 500;

% Relate the clustering performance
% Better set this value to be 200
option.zoptions.gamma = 2e10;
option.zoptions.tof = 1e-4;
option.zoptions.dis = false;
option.zoptions.distance = 'nnls';

% Initialization of variables
W = {};
H = {};
tstats = {};
indClust_NBS = {};
% Normalization of the original dataset
if flag_normalization
    X = NormalizeFea(sub_gene_data',0);
else
    X = sub_gene_data';
end
    
tic
for i = 1:n_consensus
    i
    [W{i},H{i},tstats{i}] =  nbs_nmf_cluster_network_nmf(X',...
        option.K,0,-1,option.zoptions.gamma,sparse(subnet_laplacian),...
        option.zoptions);
    % Cluster patients based on H
    indClust_NBS{i} = NMFCluster(H{i});
end
toc

% Hierarchical clustering based on the co-clustering matrix
tmp = size(sub_gene_data,2);
mat_co_nbs = zeros(tmp);
for i = 1:n_consensus
    for j = 1:tmp
        for k = 1:tmp
            if indClust_NBS{1,i}(j) == indClust_NBS{1,i}(k)
                mat_co_nbs(j,k) = mat_co_nbs(j,k)+1;
            end
        end
    end
end

cnum = option.K;
Z_nbs = linkage(squareform(n_consensus-mat_co_nbs),'av');
label_NBS = cluster(Z_nbs,'maxclust',cnum);
csvwrite('tmp.csv',mat_co_nbs);
