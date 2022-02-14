function [termWeights] = cal_term(sup_t1,Ic,lambda,gmm_U,gmm_B)
nbr_sp  = max(sup_t1(:));
termWeights = zeros(nbr_sp,2);
for i = 1:nbr_sp
termWeights_sp(i) = compute_unary(Ic(i,:), gmm_U);
termWeights_tp(i) = compute_unary(Ic(i,:), gmm_B);
end
termWeights(:,1) = lambda*real(termWeights_sp);
termWeights(:,2) = lambda*real(termWeights_tp);