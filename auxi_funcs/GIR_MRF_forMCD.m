function [ChangeMap,DifferenceImage,RegImg] = GIR_MRF_forMCD(image_t1,image_t2,opt)
%% preprocessing
Compactness = 1;
t_p1 = clock;
Compactness = 1;
[sup_img,Ns] =  SuperpixelSegmentation(image_t1,opt.Ns,Compactness);
[t1_feature,t2_feature,norm_par] = MSMfeature_extraction(sup_img,image_t1,image_t2) ;% MVE;MSM
fprintf('\n');fprintf('The computational time of Preprocessing (t_p1) is %i \n', etime(clock, t_p1));
fprintf(['\n' '====================================================================== ' '\n'])
%% Algorithm 1. Structured graph learning
opt.K = round(sqrt(Ns));
opt.iterSGL = 5;
opt.mu1 = 0.4;
opt.mu2 = 0.4;
t_p2 = clock;
[Sx] = StructuredGraphLearning(t1_feature,opt);
if strcmp(opt.Ltype,'GL')
  [GLx] = LaplacianMatrix(Sx); % GL
  elseif strcmp(opt.Ltype,'HGL')
  [GLx] = HyperLaplacianMatrix(Sx,t1_feature); % HGL
end
fprintf('\n');fprintf('The computational time of Structured graph learning (t_p2) is %i \n', etime(clock, t_p2));
fprintf(['\n' '====================================================================== ' '\n'])
%%  Algorithm 2. Structure consistency based image regression
% image_t1 ----> image_t2
opt.iterIR = 10;
opt.mu1 = 0.4;
opt.mu2 = 0.4;
t_p3 = clock;
if opt.Ns <=15000
  [regression_t1, delt,RelDiff] = StructureConsistencyImageRegression(t2_feature,GLx,Sx,opt);
elseif opt.Ns >15000 % for large Ns, the preconditioned conjugate gradient method is recommended
  [regression_t1, delt,RelDiff] = StructureConsistencyImageRegressionPCG(t2_feature,GLx,Sx,opt);
end
fprintf('\n');fprintf('The computational time of Image Regression (t_p3) is %i \n', etime(clock, t_p3));
fprintf(['\n' '====================================================================== ' '\n'])
DI_tmep = sum(delt.^2,1);
DifferenceImage  = suplabel2DI(sup_img,DI_tmep);
[RegImg,~,~] = suplabel2ImFeature(sup_img,regression_t1,size(image_t2,3));% t1--->t2
RegImg = DenormImage(RegImg,norm_par(size(image_t1,3)+1:end));
%% Algorithm 3: GMM based MRF segmentation
t_p4 = clock;
no_gauss = min(size(delt,1),fix(size(delt,2)/1000));
opt.iterGC = 5; 
% Partially use the MATLAB Implementation of GrabCut download from https://github.com/xiumingzhang/grabcut.
% And use the graphCutMex download from https://github.com/aosokin/graphCutMex_BoykovKolmogorov. 
[ChangeMap] = MeanGMMbasedMRF(image_t1,sup_img,t1_feature,t2_feature,delt,no_gauss,opt);
fprintf('\n');fprintf('The computational time of MRF segmentation (t_p4) is %i \n', etime(clock, t_p4));
fprintf(['\n' '====================================================================== ' '\n'])
