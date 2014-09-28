% This script validates COPD/TESRA 1439/3096 expression dataset using the basis matrix
% obtained from ECLIPSE 1439/3096 expression dataset

filename_W_ecl_nmf = './result_ecl/ECL_3000_Normalized_NMF_LF4/W_basis_matrix.csv';
filename_W_ecl_nbs = './result_ecl/ECL_1500_Normalized_NBS_LF4/W_basis_matrix.csv';

% NMF basis from ECL
W_ecl_nmf = csvread(filename_W_ecl_nmf);
% NBS basis from ECL
W_ecl_nbs = csvread(filename_W_ecl_nbs);

% Directory of NBS code
library_path = '/Users/changyale/matlab/nbs';
addpath(genpath(library_path))
% Add COPD dataset path
addpath('/Users/changyale/python/nbs/');
load data_gene_expression_tesra
data = NormalizeFea(data_exp_tesra,1);

% use NMF basis of ECL
H_nmf_pred = inv(W_ecl_nmf'*W_ecl_nmf)*W_ecl_nmf'*data;
label_nmf_pred = NMFCluster(H_nmf_pred);

% use NBS basis of ECL
%H_nbs_pred = inv(W_ecl_nbs'*W_ecl_nbs)*W_ecl_nbs'*data;
%label_nbs_pred = NMFCluster(H_nbs_pred);

% write results into csv file
csvwrite('H_tesra_nmf_pred.csv',H_nmf_pred);
csvwrite('label_tesra_nmf_pred.csv',label_nmf_pred);
%csvwrite('H_tesra_nbs_pred.csv',H_nbs_pred);
%csvwrite('label_tesra_nbs_pred.csv',label_nbs_pred);

