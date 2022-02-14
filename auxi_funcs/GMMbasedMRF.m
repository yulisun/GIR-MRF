function [bi_map_GC,labels] = GMMbasedMRF(sup_t1,lambda,t1_feature,t2_feature,delt_x,no_gauss,iter_max)
Ic = (sum(delt_x.^2,1));
Ic = Ic';
Ic = remove_outlier(Ic);
Ic = Ic/max(Ic);
T_theory = graythresh(Ic);
Ic_thold = (sign(Ic -T_theory) +1)/2;% 1 for foreground,0 for background
pix_U = Ic_thold==1;
Ic_feature = delt_x';
T_U = Ic_feature(pix_U, :);
pix_B = ~pix_U;
T_B = Ic_feature(pix_B, :);

addpath('GMM');
k_B = kmeans(T_B, no_gauss);
gmm_B = fit_gmm(T_B, k_B);
% Foreground
k_U = kmeans(T_U, no_gauss);
gmm_U = fit_gmm(T_U, k_U);

iter = 1;
addpath('GC');
isConverged = 0;
label_change_thold = 0.0001;
[edgeWeights] = cal_edge(sup_t1,lambda,t1_feature,t2_feature);

while iter<=iter_max && isConverged~=1
    
    
    [k_U, k_B] = assign_gauss(Ic_feature, pix_U, gmm_U, pix_B, gmm_B);
    
    
    
    [gmm_U, gmm_B] = update_gmm(Ic_feature, pix_U, k_U, pix_B, k_B);
    
    
    
    [termWeights] =  cal_term(sup_t1,Ic_feature,lambda,gmm_U,gmm_B);
    [cut, labels] = graphCutMex(termWeights, edgeWeights);
    labels_change_ratio = sum(abs(labels-double(pix_U)))/sum(labels);
    pix_U = labels==1;
    pix_B = ~pix_U;
    
    
    if labels_change_ratio <label_change_thold
        isConverged = 1;
    end
    
    idx_t1 = label2idx(sup_t1);
    for i = 1:size(t1_feature,2)
        index_vector = idx_t1{i};
        bi_map(index_vector) = labels(i);
    end
    bi_map_GC =reshape(bi_map,[size(sup_t1,1) size(sup_t1,2)]);
    iter = iter+1;
end
