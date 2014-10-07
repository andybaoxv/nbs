% This script detects whether there are outliers in three datasets:
% ECLIPSE, COPDGene, TESRA

%% Parameter Settings
flag_all = 3;
flag_normalization = true;

%%
% Load three datasets
addpath('~/python/nbs_py/');
addpath('~/dataset/ECL_CG_TESRA/')
% 1439 x 232
load data_gene_expression_ecl_1439
% 1439 x 136
load data_copd_overlap_with_ecl
% 1439 x 247
load data_tesra_overlap_with_ecl
% 3096 x 615
load data_ecl_cg_tesra
% 1439 x 615
load data_ecl_cg_tesra_1439

if flag_normalization
    data_ecl = NormalizeFea(data_gene_expression_ecl_1439,1);
    data_copd = NormalizeFea(data_copd_overlap_with_ecl,1);
    data_tesra = NormalizeFea(data_tesra_overlap_with_ecl,1);
else
    data_ecl = data_gene_expression_ecl_1439;
    data_copd = data_copd_overlap_with_ecl;
    data_tesra = data_tesra_overlap_with_ecl;
end
clear data_gene_expression_ecl_1439;
clear data_copd_overlap_with_ecl;
clear data_tesra_overlap_with_ecl;

% Apply PCA on the datasets
[coeff_ecl,score_ecl,latent_ecl] = pca(data_ecl');
[coeff_copd,score_copd,latent_copd] = pca(data_copd');
[coeff_tesra,score_tesra,latent_tesra] = pca(data_tesra');

% Apply PCA on three datasets
if flag_all == 1
    [coeff_all,score_all,latent_all] = pca([data_ecl,data_copd,data_tesra]');
elseif flag_all == 2
    [coeff_all,score_all,latent_all] = pca(data_ecl_cg_tesra');
elseif flag_all == 3
    [coeff_all,score_all,latent_all] = pca(data_ecl_cg_tesra_1439');
end


n_samples_ecl = 229;
n_samples_cg = 136;
n_samples_tesra = 247;

% Visualization on 2 dimensional space
figure(1);
% hold on
% scatter(score_ecl(:,1),score_ecl(:,2),[],[1,0,0]);
% scatter(score_copd(:,1),score_copd(:,2),[],[0,1,0]);
% scatter(score_tesra(:,1),score_tesra(:,2),[],[0,0,1]);
c = [repmat([1,0,0],[n_samples_ecl,1]);...
     repmat([0,1,0],[n_samples_cg,1]);...
     repmat([0,0,1],[n_samples_tesra,1])];
scatter(score_all(:,1),score_all(:,2),[],c);
hold off

% Visualization on 3 dimensional space
figure(2);
x = [score_ecl(:,1);score_copd(:,1);score_tesra(:,1)];
y = [score_ecl(:,2);score_copd(:,2);score_tesra(:,2)];
z = [score_ecl(:,3);score_copd(:,3);score_tesra(:,3)];
c = [repmat([1,0,0],[n_samples_ecl,1]);...
     repmat([0,1,0],[n_samples_cg,1]);...
     repmat([0,0,1],[n_samples_tesra,1])];
scatter3(score_all(:,1),score_all(:,2),score_all(:,3),[],c);
% scatter3(score_ecl(:,1),score_ecl(:,2),score_ecl(:,3),[],[1,0,0]);
% scatter3(score_copd(:,1),score_copd(:,2),score_copd(:,3),[],[0,1,0]);
% scatter3(score_tesra(:,1),score_tesra(:,2),score_tesra(:,3),[],[0,0,1]);




