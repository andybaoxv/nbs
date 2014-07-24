% This script run consensus clustering on NMF

%% Load COPDGene Expression Dataset

addpath('/home/changyale/python/nbs/');
load data_gene_expression_copd
load data_gene_expression_copd_1500
load data_gene_expression_copd_3000
load data_gene_expression_ecl
load data_gene_expression_ecl_1439

% choose which dataset to use
flag_dataset = 5;
if flag_dataset == 1
    data_use = data_exp;
elseif flag_dataset == 2
    data_use = data_exp_1500;
elseif flag_dataset == 3
    data_use = data_exp_3000;
elseif flag_dataset == 4
    data_use = data_exp_ecl;
elseif flag_dataset == 5
    data_use = data_gene_expression_ecl_1439;
else
    disp(['Error: Should Specify input dataset'])
end

% Specify the number of latent factors(#metagenes). Note that the number of
% clusters is usually set to be equal to this value in subsequent analysis
cnum = 4;

% whether to normalize the dataset
flag_normalization = true;
if flag_normalization
    X = NormalizeFea(data_use,1);
else
    X = data_use;
end

%% Add NBS code
% Directory of NBS code
library_path = '/home/changyale/matlab/nbs';
addpath(genpath(library_path))

% Number of NMF runs for consensus clustering
n_consensus = 100;

%% Parameter Settings for Standard NMF
% Maximal Number of Iterations
cczoptions.iter = 1000;
% Tolerance/Precision
cczoptions.tof = 1e-3;
% Whether to display result
cczoptions.dis = false;
% 'nnls' is an optimization technique to solve the NMF problem
cczoptions.dist = 'nnls';

%% Parameter Settings for Sparse NMF
optionnnls.eta = (max(max(X)))^2;
optionnnls.beta = 0;
optionnnls.iter = 1000;
optionnnls.dis = false;

% Initialization of variables
A = {};
Y = {};
numIter = {};
tElapsed = {};
finalResidual = {};
indClust_NMF = {};

tic
parfor i = 1:n_consensus
    % Apply NMF
    disp(['run-> ' num2str(i)]);
    %[A{i},Y{i},numIter{i},tElapsed{i},finalResidual{i}] = nmfnnls(X,cnum,cczoptions);
    [A{i},Y{i},numIter{i},tElapsed{i},finalResidual{i}] = sparsenmfnnls(X,cnum,optionnnls);
    % Normalize the columns of A{i}, make it sum to 1
    tmp_1 = diag(1./sum(A{i}));
    tmp_2 = diag(sum(A{i}));
    A{i} = A{i}*tmp_1;
    Y{i} = tmp_2*Y{i};
    % Clustering patients based on Y
    indClust_NMF{i} = NMFCluster(Y{i});  
end
toc

% Hierarchical clustering based on the co-clustering matrix
tmp = size(X,2);
mat_co_nmf = zeros(tmp);
for i = 1:n_consensus
    for j = 1:tmp
        for k = 1:tmp
            if indClust_NMF{1,i}(j) == indClust_NMF{1,i}(k)
                mat_co_nmf(j,k) = mat_co_nmf(j,k)+1;
            end
        end
    end
end

Z_nmf = linkage(squareform(n_consensus-mat_co_nmf),'av');
label_NMF = cluster(Z_nmf,'maxclust',cnum);
csvwrite('tmp.csv',mat_co_nmf);
disp(['dispersion coefficient: ' num2str(dispersionCoefficient(mat_co_nmf/n_consensus))]);
