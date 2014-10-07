% This function compute characteristic of each cluster according to
% phenotype data and labels

function [cluster_size,n_cases_cluster,n_controls_cluster,...
    fev1_mean_cluster,fev1_median_cluster,fev1_std_cluster] = ...
    compute_fev1_case_status_per_cluster(case_status,fev1,label_pred)
% Parameters
% ----------
% case_status: vector
%           case/control status for each patient
% fev1: vector
%           fev1 values for each patient
% label_pred: vector
%           predicted labels
%
% Returns
% -------

% number of clusters
label_pred_unique = unique(label_pred);
n_clusters = length(unique(label_pred));
% size of each cluster
cluster_size = zeros(n_clusters,1);
% # cases in each cluster
n_cases_cluster = zeros(n_clusters,1);
% # controls in each cluster
n_controls_cluster = zeros(n_clusters,1);
% mean of fev1 values in each cluster
fev1_mean_cluster = zeros(n_clusters,1);
% median of fev1 values in each cluster
fev1_median_cluster = zeros(n_clusters,1);
% standard deviation of fev1 values in each cluster
fev1_std_cluster = zeros(n_clusters,1);

for k = 1:n_clusters
    idx_sel = find(label_pred == label_pred_unique(k));
    cluster_size(k) = length(idx_sel);
    n_cases_cluster(k) = sum(case_status(idx_sel));
    n_controls_cluster(k) = cluster_size(k) - n_cases_cluster(k);
    fev1_mean_cluster(k) = mean(fev1(idx_sel));
    fev1_median_cluster(k) = median(fev1(idx_sel));
    fev1_std_cluster(k) = std(fev1(idx_sel));
end

% visualization of case/control status and fev1 values
figure(1);
fig_1 = bar([n_cases_cluster,n_controls_cluster],'stacked');
legend(fig_1,{'Cases','Controls'},'Location','Best','FontSize',8);
xlabel('Cluster ID');
ylabel('Number of Cases/Controls');
title('Case/Control Distribution in each cluster');

tmp = [(1:length(cluster_size))',cluster_size,n_cases_cluster,...
    fev1_mean_cluster,fev1_median_cluster,fev1_std_cluster];
csvwrite('tmp.csv',tmp);
% disp(fev1_mean_cluster');
% disp(fev1_median_cluster');
% disp(fev1_std_cluster');
% figure(2);
% bar(fev1_mean_cluster');
% xlabel('Cluster ID');
% ylabel('FEV1 Mean');
% title('FEV1 Mean in Each Cluster');
% 
% figure(3);
% bar(fev1_median_cluster');
% xlabel('Cluster ID');
% ylabel('FEV1 Median');
% title('FEV1 Median in Each Cluster');
% 
% figure(4);
% bar(fev1_std_cluster');
% xlabel('Cluster ID');
% ylabel('FEV1 Standard Deviation');
% title('FEV1 Standard Deviation in Each Cluster');

end