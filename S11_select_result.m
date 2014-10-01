% This script choose one result from multiple runs

rand_list = zeros(1,n_consensus);
for i = 1:n_consensus
    rand_list(i) = adjrand(NMFCluster(Y{i}),label_pred);
end
idx_sel = find(rand_list == max(rand_list))

disp(['Max Randlist: ' num2str(rand_list(idx_sel))]);

csvwrite('W_basis_matrix.csv',A{idx_sel(1)});
csvwrite('H_coef_matrix.csv',Y{idx_sel(1)});
csvwrite('label_pred.csv',label_pred);