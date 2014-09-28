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
% Note this is just the initial value of gamma
option.zoptions.gamma = 200;

% This parameter is very critical for the performance of the algorithm
option.zoptions.optGamma = 1;

option.zoptions.tof = 1e-4;
option.zoptions.dis = false;
option.zoptions.distance = 'nnls';

% Initialization of variables
W = {};
H = {};
tstats = {};
indClust_NBS = {};
% Normalization of the original dataset
data_use = data_gene_expression_ecl_1439;
if flag_normalization
    X = NormalizeFea(data_use,1);
else
    X = data_use;
end

load subnet_laplacian

tic
parfor i = 1:n_consensus
    disp(['run-> ' num2str(i)]);
    [W{i},H{i},tstats{i}] =  nbs_nmf_cluster_network_nmf(X,...
        option.K,0,-1,option.zoptions.gamma,sparse(subnet_laplacian),...
        option.zoptions);
    % Cluster patients based on H
    indClust_NBS{i} = NMFCluster(H{i});
end
toc

% Hierarchical clustering based on the co-clustering matrix
tmp = size(X,2);
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
