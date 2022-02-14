%%  Structured graph based image regression for unsupervised multimodal change detection
%{
Code: GIR-MRF - 2021
This is a test program for the graph based image regression and MRF segmentation method (GIR-MRF) for multimodal change detection problem.

If you use this code for your research, please cite our paper. Thank you!

Sun, Yuli, et al. "Structured graph based image regression for unsupervised multimodal change detection." 
ISPRS Journal of Photogrammetry and Remote Sensing 185 (2022): 16-31.
===================================================
%}
clear;
close all
addpath('auxi_funcs')
%% load dataset
addpath('datasets')
% Please note that the forward and backward detection results are not the same. 
% When the forward result is not satisfactory, try swapping the input order of image_t1 and image_t2 to get the backward change detection result. 
% In the future we will consider fusing the forward and backward results to improve detection performance.
dataset = '#5-Shuguang';% #1-Italy, #2-Texas, #3-Img7, #4-Img17, #5-Shuguang, #6-California, #7-Img5 
Load_dataset
fprintf(['\n Data loading is completed...... ' '\n'])
%% Parameter setting
% With different parameter settings, the results will be a little different
% Ns: the number of superpxiels,  A larger Ns will improve the detection granularity, but also increase the running time. 5000 <= Ns <= 20000 is recommended.
% lamda: sparse regularization parameter, which should be selected according to the proportion of the changed component. 
% eta: balance parameter. The smaller the lambda, the smoother the CM. 0.025<= alfa <=0.1 is recommended.
% Nseg: Number of repetitions of the MRF segmentation, averaged Nseg times.
% When testing the running time, set opt.Nseg = 1.
opt.Ns = 10000;
opt.PenaltyType= 'Fro'; % 'Fro', 'L1','L21'
opt.beta = 1;
opt.lambda = 0.1;
opt.eta = 0.025; % for #1-Italy, eta = 0.05 is better
opt.Ltype = 'GL';% 'GL' or 'HGL'
opt.Nseg = 10; % set opt.Nseg = 1 for testing running time.
%% GIR-MRF
fprintf(['\n GIR-MRF is running...... ' '\n'])
t_o = clock;
[ChangeMap,DifferenceImage,RegImg] = GIR_MRF_forMCD(image_t1,image_t2,opt);

fprintf('\n');fprintf('The total computational time of GIR-MRF (t_total) is %i \n', etime(clock, t_o));
%% Displaying results
fprintf(['\n' '====================================================================== ' '\n'])
fprintf(['\n Displaying the results...... ' '\n'])
%---------------------AUC PCC F1 KC ----------------------%
n=500;
Ref_gt = Ref_gt/max(Ref_gt(:));
[TPR, FPR]= Roc_plot(DifferenceImage,Ref_gt,n);
[AUC, Ddist] = AUC_Diagdistance(TPR, FPR);
for i = 1:size(ChangeMap,3)
    [tp,fp,tn,fn,fplv,fnlv,~,~,pcc(i),kappa(i),imw]=performance(ChangeMap(:,:,i),Ref_gt);
    F1(i) = 2*tp/(2*tp + fp + fn);
end
result = 'AUC is %4.3f; PCC is %4.3f; F1 is %4.3f; KC is %4.3f \n';
[~,idx_best] = max(F1);
fprintf(result,AUC,pcc(idx_best),F1(idx_best),kappa(idx_best))

%------------Regression image,  Difference imag and Change map --------------%
figure; plot(FPR,TPR);title('ROC curves');
figure;
if strcmp(dataset,'#2-Texas') == 1
    subplot(131);imshow(6*uint16(RegImg(:,:,[5 4 3])));title('Regression image')
elseif strcmp(dataset,'#6-California') == 1
    subplot(131);imshow(RegImg,[-1 1]);title('Regression image')
elseif strcmp(dataset,'#7-Img5') == 1
    subplot(131);imshow(uint8(exp(RegImg*3.75+1.8)));title('Regression image')
else
    subplot(131);imshow(uint8(RegImg));title('Regression image')
end
subplot(132);imshow(remove_outlier(DifferenceImage),[]);title('Difference image')
subplot(133);imshow(ChangeMap(:,:,idx_best),[]);title('Change mape')
