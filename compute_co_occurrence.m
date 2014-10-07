% This function compute co_occurrence matrix from clustering solutions of
% multiple runs

function mat_co = compute_co_occurrence(ind_cluster)

n_consensus = length(ind_cluster);
n_samples = length(ind_cluster{1,1});
mat_co = zeros(n_samples);

for i = 1:n_consensus
    for j = 1:n_samples
        for k = 1:n_samples
            if ind_cluster{1,i}(j) == ind_cluster{1,i}(k)
                mat_co(j,k) = mat_co(j,k)+1;
            end
        end
    end
end

end