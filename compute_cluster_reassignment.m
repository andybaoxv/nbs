% This script show the change of cluster assignment when the number of
% latent factors changes from 3 to 4

load label_nbs_lf3
% Order from LF3 NBS
load Order
load label_nbs_lf4
load label_nbs_lf5
load label_nbs_cg_lf3
load label_nbs_cg_lf4
load label_nbs_cg_lf5
load label_nbs_tesra_lf3
load label_nbs_tesra_lf4
load label_nbs_tesra_lf5

label_nbs_lf3_order = label_nbs_lf3(Order);
label_nbs_lf4_order = label_nbs_lf4(Order);
label_nbs_lf5_order = label_nbs_lf5(Order);

res_ecl = [label_nbs_lf3_order,label_nbs_lf4_order,label_nbs_lf5_order];
figure(1)
imagesc(res_ecl);
colorbar

res_cg = [label_nbs_cg_lf3(Order),label_nbs_cg_lf4(Order),...
    label_nbs_cg_lf5(Order)];
figure(2)
imagesc(res_cg);
colorbar

res_tesra = [];
