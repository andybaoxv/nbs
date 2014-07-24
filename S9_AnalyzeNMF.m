% This script analyze clustering solution from NMF on gene expression data
flag_normalization_pheno = true;

W = csvread('NMF_norm_LF4_basis.csv');
H = csvread('NMF_norm_LF4_coef.csv');
[W_sorted,idx_W] = sort(W,'descend');


label_pred = NMFCluster(H);
label_pred_unique = unique(label_pred);
n_features = size(data,2);
n_clusters = length(label_pred_unique);

load data_pheno
load patient_label
if flag_normalization_pheno
    data = normalization(data_pheno);
else
    data = data_pheno;
end
    
cluster_data = {};

% Separate original data by cluster
for i = 1:n_clusters
    idx_tmp = find(label_pred == label_pred_unique(i));
    cluster_data{label_pred_unique(i)} = data(idx_tmp,:);
end

%% Compute statistics for each cluster
% mean of each cluster
cluster_mean = zeros(n_clusters,n_features);
cluster_std = zeros(n_clusters,n_features);
cluster_min = zeros(n_clusters,n_features);
cluster_max = zeros(n_clusters,n_features);
cluster_med = zeros(n_clusters,n_features);
cluster_case = zeros(n_clusters,1);
cluster_control = zeros(n_clusters,1);

for k = 1:n_clusters
    idx_tmp = find(label_pred == label_pred_unique(k));
    cluster_mean(k,:) = mean(cluster_data{k});
    cluster_std(k,:) = std(cluster_data{k});
    cluster_min(k,:) = min(cluster_data{k});
    cluster_max(k,:) = max(cluster_data{k});
    cluster_med(k,:) = median(cluster_data{k});
    cluster_case(k) = sum(patient_label(idx_tmp)==1);
    cluster_control(k) = sum(patient_label(idx_tmp)==0);
end

if flag_normalization_pheno
    csvwrite('cluster_mean_norm.csv',cluster_mean);
    csvwrite('cluster_std_norm.csv',cluster_std);
    csvwrite('cluster_min_norm.csv',cluster_min);
    csvwrite('cluster_max_norm.csv',cluster_max);
    csvwrite('cluster_med_norm.csv',cluster_med);
else
    csvwrite('cluster_mean_orig.csv',cluster_mean);
    csvwrite('cluster_std_orig.csv',cluster_std);
    csvwrite('cluster_min_orig.csv',cluster_min);
    csvwrite('cluster_max_orig.csv',cluster_max);
    csvwrite('cluster_med_orig.csv',cluster_med);
end
csvwrite('cluster_case.csv',cluster_case);
csvwrite('cluster_control.csv',cluster_control);
% tmp = [cluster_case,cluster_control];
% figure
% bar(tmp,'stacked');
% xlabel('Cluster ID');
% ylabel('Number of Case/Control');
% title('Case(blue)/Control(red) Distribution in Each Cluster');

% features_name = {'BD.FEV1','oxygen','ExacTrunc','BD.FEV.FVC',...
%     'FracVol.950U','Lowest.15.','Emphysema','Neutrophils','Lymphocytes',...
%     'Monocytes','Eosinophils','Basophils'};
% 
% figure
% tmp = bar3(cluster_med);
% hx = get(tmp(1),'parent');
% set(hx,'xtickLabel',features_name);
% set(hx,'ytickLabel',{'Clust 1','Clust 2','Clust 3','Clust 4'});
