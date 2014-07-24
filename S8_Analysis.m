% Compute NMI between different solutions

load label_corr_norm
load label_corr_pca
load label_nmf_lf4
load label_nbs_lf4
load label_nbs_lf2
load label_nbs_lf3
load label_nbs_lf10
load patient_label
load data_pheno

for i = 1:length(patient_label)
    if patient_label(i) == 0
        patient_label(i) = 2;
    end
end

sz_pheno = size(data_pheno);
for j = 1:length(sz_pheno(2))
    nmi(label_nbs_lf4,data_pheno(:,j))
end

