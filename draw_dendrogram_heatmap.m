% This function draw heatmap and denrogram given similarity matrix
function [fig_1,fig_2] = draw_dendrogram_heatmap(S,idx_fig_1,idx_fig2)
% Parameters
% ----------
% S: matrix, shape(n_samples,n_samples)
%    similarity matrix
%
% idx_fig_1,idx_fig_2: integers
%    indices of figure id
%
% Returns
% --------
% fig_1,fig_2: handlers for generated figures

% number of samples
[num_sample,~] = size(S);

% normalize similarity matrix
S_max = max(max(S));
S = S/S_max;

% B is the corresponding distance matrix
B = 1-S;
Dist = [];
for m = 1:(num_sample-1)
    Dist = cat(2,Dist,B(m,(m+1):num_sample));
end
orig_tree  = seqlinkage(Dist,'average');

% draw dendrogram with linkage
Tree = linkage(Dist,'average');
fig_1 = figure(idx_fig_1);
% Order = optimalleaforder(Tree, Dist);
[H_1,T,perm] = dendrogram(Tree,0,'orientation','top');

% draw heatmap
Order = perm;
C = zeros(num_sample,num_sample);
for m = 1:num_sample
    for n = 1:num_sample
        C(m,n) = S(Order(m),Order(n));
    end
end

fig_2 = figure(idx_fig_2);
imagesc(C,[0,1]);

end