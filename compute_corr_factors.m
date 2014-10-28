% This script compute pairwise corr between factors in basis matrix

mtr_factor = W_sel;

n_factors = size(mtr_factor,2);
corr_factor = zeros(n_factors,n_factors);

for i = 1:n_factors
    for j = i:n_factors
        corr_factor(i,j) = corr(mtr_factor(:,i),mtr_factor(:,j));
        corr_factor(j,i) = corr_factor(i,j);
    end
end

disp(corr_factor);