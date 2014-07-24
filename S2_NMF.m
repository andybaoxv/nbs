% This script apply NMF on the extracted dataset

%% Load COPDGene Expression Dataset
addpath('/Users/changyale/python/nbs/');
load data_gene_expression_copd
sub_gene_data = data_exp;
clear data_exp;

%% Parameter Settings for Standard NMF
% whether to normalize dataset
flag_normalization = true;
% number of latent factors
cnum = 4;
% parameters for NMF
optionsnmf.alpha2 = 1e-5;
optionsnmf.alpha1 = 1e-5;
optionsnmf.lambda2 = 0;
optionsnmf.lambda1 = 0;
optionsnmf.t1 = true;
optionsnmf.t2 = true;
optionsnmf.iter=1000;
optionsnmf.dis=false;
optionsnmf.residual=1e-4;
optionsnmf.tof=1e-4;

% Data normalization
if flag_normalization
    %X = normc(sub_gene_data);
    X = NormalizeFea(sub_gene_data,1);
else
    X = sub_gene_data;
end

% Sparse NMF Settings
optionnnls.eta = (max(max(X)))^2;
optionnnls.beta = 0;
optionnnls.iter = 1000;
optionnnls.dis = false;

% NMF with regularization
addpath('/Applications/MATLAB_R2014a.app/toolbox/nmfv1_4');
A_list = {};
Y_list = {};
num_run = 10;
finalResidual_list = zeros(num_run,1);
for i = 1:num_run
    disp(['run-> ' num2str(i)]);
    %[A,Y,AtA,numIter,tElapsed,finalResidual_list(i)]=vsmf(X,cnum,optionsnmf);
    [A,Y,numIter,tElapsed,finalResidual_list(i)]=sparsenmfnnls(X,cnum,optionnnls);
    A_list{i} = A;
    Y_list{i} = Y;
end
idx_min = find(finalResidual_list == min(finalResidual_list));
idx_sel = idx_min(1);
A = A_list{idx_sel};
Y = Y_list{idx_sel};

%A = A_sel;
%Y = Y_sel;
%label = label_sel;

% gene ranking by gene score
optionffnmf.isBasis=true;
optionffnmf.propertyName='median';
optionffnmf.propertyValue=1;
[mask1,featNames,scores1,A1]=featureFilterNMF(A,[],optionffnmf);

% for each basis vector, find the factor-specific genes


for i=1:cnum
    basisVec=A(:,i);
    maxVal=max(basisVec);
    maski=(basisVec>1e-3*maxVal); % nonzero logical indicators
    mask_max = (basisVec == max(A')');
    tmp = mask1 & maski & mask_max;
    genei=sub_gene_hgnc_symbol(tmp); % genes for the i-th factor
    fprintf('The number of genes in the %d-th factor is %d\n.',i,numel(genei));
    fileName=['./Factor',num2str(i)];
    writeDir='factors';
    writeGeneList(genei,fileName,writeDir); % write into txt file
end


% %% Standard NMF
% 
% % INPUT: gene_data: p x n data matrix;(p=31, n=121 in our case)
% %        cnum: the number of latent factors(megagenes)
% %        cczoptions: parameter settings
% % OUTPUT: A: p x q basis vector matrix
% %         Y: q x n coefficient matrix(low dimensional embedding matrix)
% cnum_range = 2:50;
% % Rank of the matrix consisting of basis vectors
% rank_A = zeros(1,length(cnum_range));
% nuclear_A = zeros(1,length(cnum_range));
% for i = 1:length(cnum_range)
%     cnum = cnum_range(i)
%     [A,Y,numIter,tElapsed,finalResidual] = nmfnnls(X',cnum,cczoptions);
%     rank_A(i) = rank(A);
%     nuclear_A(i) = sum(svd(A));
% end
% subplot(2,1,1);
% plot(cnum_range,rank_A);
% subplot(2,1,2);
% plot(cnum_range,nuclear_A);
% % %% Clustering patients based on obtained low dimensional matrix Y
% 
% INPUT: Y: q x n coefficient matrix. Note that each column is a sample.
% OUTPUT: indClust_ctrl: n x 1 cluster membership indicator

%indClust_ctrl = NMFCluster(Y);
indClust_ctrl = kmeans(Y',cnum);
counter = zeros(1,cnum);
for i = 1:length(indClust_ctrl)
    counter(indClust_ctrl(i)) = counter(indClust_ctrl(i))+1; 
end
counter

