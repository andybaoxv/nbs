clear all, close all, clc;

tic
%% Load STRING gene interaction network contained in NBS data folder

% Directory of NBS data
basedata_path = '/Users/changyale/dataset/nbs_data';

% Directory of NBS code
library_path = '/Users/changyale/matlab/nbs';
addpath(genpath(library_path))

% load STRING network, which contain 12233 genes.
NBSnetwork = load([basedata_path '/networks/ST90Q_adj_mat.mat']);

%% Load an NBS influence matrix

% Influence matrix is derived from the original STRING network. The
% motivation is to build a more accurate measure of similarity between
% genes. The size of this matrix is 12233 x 12233
glap = load([basedata_path '/networks/glap_subnetwork_ST90.mat']);

%% Load Eclipse Gene Expression Data 

%load gene_eclipse.mat
load eclipse3.mat
load patient_label.mat

gene_eclipse = eclipse3;

[n_row,n_col] = size(gene_eclipse);

% Extract Gene Expression Data
gene_data = cell2mat(gene_eclipse(2:n_row,5:n_col));

% Extract patient ID
tmp = gene_eclipse';
patient_id = tmp(5:n_col,1);

% Extract probe names
gene_probe_name = gene_eclipse(2:n_row,1);

% Extract HGNC Gene Symbol
gene_hgnc_symbol = gene_eclipse(2:n_row,2);

% Extract Ensembl Gene ID
gene_ensembl_id = gene_eclipse(2:n_row,3);

% Extract Entrez Gene ID
gene_entrez_id = cell2mat(gene_eclipse(2:n_row,4));

% tmp
tmp = {zeros(n_col-4,1)};
tmp_1 = struct('center',tmp,'tss',tmp,'filesource',tmp);
% Construct a struc similar to baseMData in demo_NBS_TCGA.m
tp_1 = struct('gene_id_all',gene_entrez_id);
tp_2 = struct('gene_indiv_mat',gene_data');
tp_3 = struct('sample_id',{patient_id});
tp_4 = struct('batchls',tmp_1);
tp_5 = struct('gene_id_symbol',{gene_hgnc_symbol});

gene_exp = struct('gene_id_all',gene_entrez_id,'gene_indiv_mat',gene_data',...
    'sample_id',{patient_id},'batchls',tmp_1,'gene_id_symbol',{gene_hgnc_symbol});
clear tp_1;
clear tp_2;
clear tp_3;
clear tp_4;
clear tp_5;
clear tmp;
clear tmp_1;


%% Construct Graph Laplacian from loaded influence matrix

% The number of nearest neighbors in constructing matrix K from graph
% Laplacian of influence matrix.
%outDeg = 11;

% glap.glap is a influence matrix, which is also a similarity matrix
% measuring pairwise similarity between genes.
% Step (1) Construct influence matrix S by defining diffusion process in
% the original gene interaction network 
% Step (2) Construct knnMat from S by setting the number of nearest 
% neighbors, which is denoted by outDeg.
% Step (3) Return the Graph Laplacian of knnMat and get knnGlap, which is
% the matrix K in NBS paper. 
% Note that the degree to which W'KW respect local gene graph topology is
% controlled by the number of nearest neighbors in constructing K, i.e.
% the value of outDeg. 
% Also note that matrix K becomes very sparse while the original influence 
% matrix is not sparse at all. (There's no 0 in glap, the percentage of 
% nonzero elements in knnGlap is 0.12%)
% It's interesting to see why use KNN instead of influence matrix directly.
%[knnGlap] = sim_to_knn_glap(glap.glap,outDeg);

toc