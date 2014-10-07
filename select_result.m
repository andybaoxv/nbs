% This function select basis matrix and coefficient matrix from multiple
% runs of NMF and NBS

function [W,H] = select_result(W_list,H_list,label_pred)

% number of consensus runs
n_consensus = length(W_list);
% rand index values for the clustering result of each run
rand_list = zeros(1,n_consensus);
for i = 1:n_consensus
    rand_list(i) = adjrand(NMFCluster(H_list{i}),label_pred);
end
idx_sel = find(rand_list == max(rand_list));

W = W_list{idx_sel(1)};
H = H_list{idx_sel(1)};

end