% This script choose one result from multiple runs

rand_list = zeros(1,n_consensus);
for i = 1:n_consensus
    rand_list(i) = adjrand(NMFCluster(Y{i}),label_NMF);
end
idx_sel = find(rand_list == max(rand_list))

