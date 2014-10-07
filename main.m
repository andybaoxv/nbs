% This script run consensus clustering on NMF and NBS
clear all, close all, clc;

%% Parameter Settings
% choose which dataset to use: 1->3096, 2->1439
flag_dataset = 2;
% whether to normalize dataset
flag_normalization = true;
% number of runs for consensus clustering
n_consensus = 100;

% ************** parameter settings for standard NMF *******************
% number of latent factors for NMF
cnum = 4;
% maximal number of iterations
options_nmf.iter = 1000;
% Tolerance/Precision
options_nmf.tof = 1e-4;
% Whether to display result
options_nmf.dis = false;
% 'nnls' is an optimization technique to solve the NMF problem
options_nmf.dist = 'nnls';
% *************** parameter settings for NBS ************************
% Set the number of latent factors for NBS
options_nbs.K = 4;
options_nbs.iter = 500;
options_nbs.zoptions.iter = 500;
% initial value of gamma
options_nbs.zoptions.gamma = 200;
% This parameter is very critical for the performance of the algorithm
options_nbs.zoptions.optGamma = 1;
options_nbs.zoptions.tof = 1e-4;
options_nbs.zoptions.dis = false;
options_nbs.zoptions.distance = 'nnls';

%% Load Gene Expression Dataset for three cohorts after COMBAT
% directory of the original dataset
addpath('~/dataset/ECL_CG_TESRA/');
% directory of NBS code
library_path = '~/matlab/nbs';
addpath(genpath(library_path))

% 3096 x 612 (ECL_229 + CG_136 + TESRA_247), #probes = 3096
load data_ecl_cg_tesra
% 1439 x 612 (ECL_229 + CG_136 + TESRA_247), #probes = 1439
load data_ecl_cg_tesra_1439
% Phenotype info for these 612 patients
load case_status_ecl_cg_tesra
load fev1_values_ecl_cg_tesra
% Load laplacian matrix for 1439 probes
load subnet_laplacian
% Number of patients in each cohort
n_ecl = 229;
n_cg = 136;
n_tesra = 247;

% use 3096 probes version
if flag_dataset == 1
    data_use = data_ecl_cg_tesra;
% use 1439 probes version
elseif flag_dataset == 2
    % separate 612 patients into three cohorts: ECL, CG, TESRA
    % gene expression data
    data_ecl_1439 = data_ecl_cg_tesra_1439(:,1:n_ecl);
    data_cg_1439 = data_ecl_cg_tesra_1439(:,(n_ecl+1):(n_ecl+n_cg));
    data_tesra_1439 = data_ecl_cg_tesra_1439(:,(n_ecl+n_cg+1):n_ecl+n_cg+n_tesra);
    % case status
    case_status_ecl = case_status_ecl_cg_tesra(1:n_ecl);
    case_status_cg = case_status_ecl_cg_tesra((n_ecl+1):(n_ecl+n_cg));
    case_status_tesra = case_status_ecl_cg_tesra((n_ecl+n_cg+1):(n_ecl+n_cg+n_tesra));
    % fev1 values
    fev1_ecl = fev1_values_ecl_cg_tesra(1:n_ecl);
    fev1_cg = fev1_values_ecl_cg_tesra((n_ecl+1):(n_ecl+n_cg));
    fev1_tesra = fev1_values_ecl_cg_tesra((n_ecl+n_cg+1):(n_ecl+n_cg+n_tesra));
    % data used for applying NMF
    data_use = data_ecl_1439;
else
    disp(['Error: Should Specify input dataset'])
end

% whether to normalize the dataset, here normalize the variance of each
% probe to be 1
if flag_normalization
    X = NormalizeFea(data_use,1);
else
    X = data_use;
end

% Initialization of variables for NMF
% A is basis matrix for NMF
A = {};
% Y is coefficient matrix for NMF
Y = {};
numIter = {};
tElapsed = {};
finalResidual = {};
indClust_NMF = {};
% initialization of variables for NBS
% W is basis matrix for NBS
W = {};
% H is coefficient matrix for NBS
H = {};
tstats = {};
indClust_NBS = {};

% use parallel computing
parfor i = 1:n_consensus
    % apply NMF
    disp(['run-> ' num2str(i)]);
    [A{i},Y{i},numIter{i},tElapsed{i},finalResidual{i}] = nmfnnls(X,cnum,options_nmf);
    % normalize the columns of A{i}, make it sum to 1
    tmp_1 = diag(1./sum(A{i}));
    tmp_2 = diag(sum(A{i}));
    A{i} = A{i}*tmp_1;
    Y{i} = tmp_2*Y{i};
    % cluster samples based on coefficient matrix Y
    indClust_NMF{i} = NMFCluster(Y{i});
    
    % apply NBS
    [W{i},H{i},tstats{i}] =  nbs_nmf_cluster_network_nmf(X,...
        options_nbs.K,0,-1,options_nbs.zoptions.gamma,sparse(subnet_laplacian),...
        options_nbs.zoptions);
    % cluster samples based on coefficient matrix H
    indClust_NBS{i} = NMFCluster(H{i});
end

%% compute co-occurrence matrix for NMF and NBS
mat_co_nmf = compute_co_occurrence(indClust_NMF);
mat_co_nbs = compute_co_occurrence(indClust_NBS);

% hierarchical clustering for NMF and NBS
Z_nmf = linkage(squareform(n_consensus-mat_co_nmf),'av');
label_nmf = cluster(Z_nmf,'maxclust',cnum);
Z_nbs = linkage(squareform(n_consensus-mat_co_nbs),'av');
label_nbs = cluster(Z_nbs,'maxclust',options_nbs.K);

% draw dendrograms and heatmaps for nmf
draw_dendrogram_heatmap(mat_co_nmf,1,2);
% draw dendrograms and heatmaps for nbs
draw_dendrogram_heatmap(mat_co_nbs,3,4);

% select result from n_consensus runs for NMF and NBS
[A_sel,Y_sel] = select_result(A,Y,label_nmf);
[W_sel,H_sel] = select_result(W,H,label_nbs);

