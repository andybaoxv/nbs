% This script will find the indices in STRING network that correspond to
% genes in the p genes in our gene expression data.

tic
%% Extract id by comparing ensembl ID

% gene_ensembl_id - ensembl id of genes in expression data, len=p
% net_ensembl_id - ensembl id of genes in the whole network, len=d=12233
net_ensembl_id = NBSnetwork.network.key;

% Overlap between genes in gene ensembl id and STRING network ensembl id 
keylist = {};

% indices of genes in expression dataset appearing in the network
keylist_gene_ind = [];

% indices of genes in the network corresponding to keylist_gene_ind
keylist_net_ind = [];

% Number of elments in common
ind = 1;

for i = 1:length(gene_ensembl_id)
    for j = 1:length(net_ensembl_id)
        if strcmp(gene_ensembl_id(i),net_ensembl_id(j))
            keylist{ind} = gene_ensembl_id{i};
            keylist_gene_ind(ind) = i;
            keylist_net_ind(ind) = j;
            ind = ind+1;
        end
    end
end

%% Extract subnetwork from STRING and graph Laplacian of influence matrix 
% through ensembl ID

% Dimension of the extracted subnetwork
dim_subnetwork = ind - 1;

% subnetwork of STRING, use NBSnetwork.network.adj_mat_norm
% After checking, this matrices contains only 0 value elments
% So it's not likely we can do smoothing anymore
subnet_string = NBSnetwork.network.adj_mat_norm(keylist_net_ind,keylist_net_ind);

% subnetwork of influence matrix
% Note that this is a symmetric matrix
% Number of nearest neighbors
n_neighbors = 11;
subnet_influence = glap.glap(keylist_net_ind,keylist_net_ind);
subnet_laplacian = sim_to_knn_glap(subnet_influence,n_neighbors);
% percentage of non-zeros elements in subnet_laplacian
%perc_nonzero = sum(sum(subnet_laplacian ~= 0))/(1439*1439);
% % Degree matrix and Lalacian of influence matrix
% subnet_degree = diag(sum(subnet_influence));
% subnet_laplacian = subnet_degree - subnet_influence;

%% Extract subset from original gene expression data

% extract subset of data
sub_gene_data = gene_data(keylist_gene_ind,:);

% extract subset of entrez id
sub_gene_entrez_id = gene_entrez_id(keylist_gene_ind);

% extract subset of HGNC symbols
sub_gene_hgnc_symbol = gene_hgnc_symbol(keylist_gene_ind);

toc