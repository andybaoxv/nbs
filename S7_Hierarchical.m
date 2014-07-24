%% This script apply hiearchical clustering on original gene expression data
% in order to cluster patients

%% Parameter Settings
% Dataset to use: 'full'(2212 x 232), 'sub'(1439 x 232)
flag_dataset = 'full';
% Whether to apply PCA before hierarchical clustering
flag_pca = false;
% Whether to select genes with top variances
flag_variance = false;
% Whether to do normalization for each feature
flag_normalization = true;
% Distance Metric to use in hierarchical clustering, including 
% 'correlation','euclidean','cosine'
flag_metric = 'correlation';

%%
% Gene expression data before mapping to STRING network
[n_gene,~] = size(gene_data);

% Gene expression data after mapping to STRING network
[n_gene_sel,~] = size(sub_gene_data);

% Choose which dataset to use, note data_use is samples x features
if strcmp(flag_dataset,'full')
    data_use = gene_data';
elseif strcmp(flag_dataset,'sub')
    data_use = sub_gene_data';
end

% Apply PCA before hierarchical clustering to avoid the curse of
% dimensionality.
if flag_pca
    [coeff,score,latent,tsquared,explained,mu] = pca(data_use);
    n_components = 0;
    for i = 1:length(latent)
        if sum(latent(1:i))/sum(latent) > 0.9
            n_components = i;
            break
        end
    end
    data_use = score(:,1:n_components);
end

% The number of genes and the number of samples
[n_instances,n_features] = size(data_use);

% Select genes with top variance
if flag_variance
    features_sel = {};
    index = 1;
    var_features = zeros(1,n_features);
    for j = 1:n_features
        var_features(j) = std(data_use(:,j));
        if var_features(j)>0.8
            features_sel{index} = j;
            index = index+1;
        end
    end
    features_sel = cell2mat(features_sel);
    data_use = data_use(:,features_sel);
    [n_instances,n_features] = size(data_use);
end

% Data Normalization, treat each gene as a variable
if flag_normalization
   for j = 1:n_features
       tmp = data_use(:,j);
       data_use(:,j) = (tmp-mean(tmp))/std(tmp);
   end
end

%% Hierarchical Clustering
% dendrogram with linkage
Tree = linkage(data_use,'average',flag_metric);
fid1 = figure(1);
[H,T,perm] = dendrogram(Tree,0,'orientation','left');
Order = perm;
% Heatmap according to hiearchical clustering result
tmp = pdist2(data_use,data_use,flag_metric);
mtr_sim = 1-tmp/max(max(tmp));
C = zeros(n_instances,n_instances);
for m = 1:n_instances
    for n = 1:n_instances
        C(m,n) = mtr_sim(Order(m),Order(n));
    end
end

label_pred = cluster(Tree,'maxclust',2);
idx_clust_1 = find(label_pred == 1);
idx_clust_2 = find(label_pred == 2);
data_clust_1 = data_use(idx_clust_1,:);
data_clust_2 = data_use(idx_clust_2,:);
% Apply PCA
[coeff_1,score_1,latent_1] = pca(data_clust_1);
for i = 1:length(latent_1)
    if sum(latent_1(1:i))/sum(latent_1) > 0.8
        num_pc_1 = i;
        break
    end
end

[coeff_2,score_2,latent_2] = pca(data_clust_2);
for i = 1:length(latent_2)
    if sum(latent_2(1:i))/sum(latent_2) > 0.8
        num_pc_2 = i;
        break
    end
end

% Compute the correlation between each gene and predicted labels
corr_gene = zeros(n_features,1);
for i = 1:n_features
    corr_gene(i) = corr(label_pred,data_use(:,i));
end
[corr_gene_sorted,idx_gene] = sort(corr_gene,'descend');

data_top_gene = data_use(:,idx_gene(1:80));
hgnc_top_gene = gene_hgnc_symbol(idx_gene(1:80));

load patient_label
label_true = patient_label+1;
nmi(label_pred,label_true)

fid2 = figure(2);
imagesc(C,[0,1]);
