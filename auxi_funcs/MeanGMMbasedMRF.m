function [CM] = MeanGMMbasedMRF(image_t1,sup_img,t1_feature,t2_feature,delt,no_gauss,opt)
[M,N,~] = size(image_t1);
CM = zeros(M,N,opt.Nseg);
error_num = 0;
for i = 1:opt.Nseg
    try
        [gc_result,~] = GMMbasedMRF(sup_img,opt.eta,t1_feature,t2_feature,delt,no_gauss,opt.iterGC);
        CM(:,:,i) = gc_result;
    catch
        error_num=error_num+1;
        continue
    end
end
idx_err= [];
ii = 1;
for i = 1:opt.Nseg
    if sum(sum(CM(:,:,i)))==0
         idx_err(ii) = i;
         ii = ii+1;
    elseif sum(sum(CM(:,:,i)))==M*N
         idx_err(ii) = i;
         ii = ii+1;
    end
end
idx_right = setdiff([1:opt.Nseg],idx_err);
CM = CM(:,:,idx_right);