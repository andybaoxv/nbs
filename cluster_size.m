% This function computes cluster size

function counter_NMF_CC = cluster_size(label_NMF)

% Counte the cluster size
counter_NMF_CC = zeros(1,4);
for i = 1:length(label_NMF)
    counter_NMF_CC(label_NMF(i)) = counter_NMF_CC(label_NMF(i))+1; 
end

end