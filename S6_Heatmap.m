
% Load similarity matrix
S = csvread('tmp.csv');
[num_sample,~] = size(S);
A_max = max(max(S));

% S is normalized similarity matrix, the values in A is between 0 an 1
S = S/A_max;

% B is the corresponding distance matrix
B = 1-S;
Dist = [];
for m = 1:(num_sample-1)
    Dist = cat(2,Dist,B(m,(m+1):num_sample));
end
orig_tree  = seqlinkage(Dist,'average');

%% dendrogram with linkage
Tree = linkage(Dist,'average');
fid1 = figure(1);
% Order = optimalleaforder(Tree, Dist);
[H_1,T,perm] = dendrogram(Tree,0,'orientation','left');
Order = perm;

C = zeros(num_sample,num_sample);
for m = 1:num_sample
    for n = 1:num_sample
        C(m,n) = S(Order(m),Order(n));
    end
end

fid2 = figure(2);
imagesc(C,[0,1]);

label_pred = cluster(Tree,'maxclust',cnum);

% tp_label_gold = zeros(n_instances-sum(gold_stage==0),1);
% tp_label_pred = zeros(n_instances-sum(gold_stage==0),1);
% idx = 1;
% for i = 1:n_instances
%     if gold_stage(i) ~= 0
%         tp_label_gold(idx) = gold_stage(i);
%         tp_label_pred(idx) = label_pred(i);
%         idx = idx+1;
%     end
% end

% set(gca,'YTick',[1:num_sample]);
% set(gca,'YTickLabel',label_verify);
% set(gca,'XTick',[]);
% colorbar;
% saveas(fid2,'heatmap.png')
%saveas(fid2,'fig_2','fig')