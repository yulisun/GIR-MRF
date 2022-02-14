function [W] = StructuredGraphLearning(X,opt)
%min{L(E,W,S)=2Tr(X*L*X')+beta*g(E)+alfa||W||_F^2} 
%s.t. W=S, E=X-XW,W1=1,0<W<1;
[M,N]=size(X);
%% adaptive K
K=opt.K;
K = K+1;
[idx,~] = knnsearch(X',X','k',K);
degree = tabulate(idx(:));
Kmat = degree(:,2);
kmax = K;
kmin = round(kmax/10)+1;
Kmat(Kmat>=kmax)=kmax;
Kmat(Kmat<=kmin)=kmin;
if length(Kmat)<N
    Kmat(length(Kmat)+1:N) = kmin;
end
%% W
S = zeros(N,N);
mu1 = opt.mu1;
mu2 = opt.mu2;
beta = opt.beta;
iter = opt.iterSGL;
R1 = zeros(N,N);
R2 = zeros(M,N);
DistX_temp = pdist(X');
DistX_temp = DistX_temp.^2;
DistX = squareform(DistX_temp);
inv_IXTX=eye(N)-X'*inv(mu1*eye(M)/mu2+X*X')*X;
for ii=1:iter
    if ii>=2
        Sold=W;
    end
    [E] = ErrorUpdate(R2,X,S,beta,mu2,opt.PenaltyType);
    P = DistX+R1-mu1*S;
    [Pasc,idxP]=sort(P);
    W = zeros(N,N);
    for i = 1:N
        K = Kmat(i);
        k = K-1;
        id_x = idxP(1:K,i);
        di = Pasc(1:K,i);
        s_w = (di(K)-di)/(k*di(K)-sum(di(1:k))+eps);
        W(id_x,i) = s_w;
    end
    S = inv_IXTX*(W+(mu2/mu1)*X'*(X-E)+(R1+X'*R2)/mu1);
    R1 = R1+mu1*(W-S);
    R2 = R2+mu2*(X-X*S-E);
    if ii>3 && norm(W-Sold,'fro')/norm(Sold,'fro')<1e-3
        break
    end
end










