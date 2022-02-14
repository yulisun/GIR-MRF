function [Z, delt,RelDiff] = StructureConsistencyImageRegression(Y,L,W,opt)
%min{L(Z,delt,E)=2Tr(ZLZ')+beta*g(E)+lamda||delt||2,1} 
%s.t. Y=Z-delt, E=Z-ZW
% X----> Y
beta = opt.beta;
iter = opt.iterIR;
lamda = opt.lambda;
mu1 = opt.mu1;
mu2 = opt.mu2;
[M,N] = size(Y);
delt = zeros(M,N);% initialization
R1 = zeros(M,N);% initialization
R2 = zeros(M,N);
IZ_IZT = mu2*((eye(N)-W)*(eye(N)-W)');
IZT = (eye(N)-W)';
inv_Tmu = inv(mu1*eye(N)+4*L+IZ_IZT);
Z = Y;

for i=1:iter
    delt_old = delt;
    [E] = ErrorUpdate(R2,Z,W,beta,mu2,opt.PenaltyType);% g(E)
    Q = Z-Y+R1/mu1;
    delt = deltUpdate(Q,lamda/mu1,21);% delt update£» 1---> L1 norm; 21---> L21 norm
    Z = (mu1*(delt+Y)+mu2*E*IZT-R1-R2*IZT)*inv_Tmu; % Z update
    R1 = R1 + mu1*(Z-Y-delt); % R1 update
    R2 = R2 + mu2*(Z-Z*W-E); % R2 update
    RelDiff(i) = norm(delt-delt_old,'fro')/norm(delt,'fro');
    if i>3 && RelDiff(i)<1e-2
        break
    end
end