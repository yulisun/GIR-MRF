function [error]= ErrorUpdate(R2,X,S,beta,mu,PenaltyType)
switch PenaltyType
    case 'L1'
        Q = R2/mu+(X-X*S);
        E = zeros(size(Q));
        epsilon = beta/mu;
        DD = abs(Q)-epsilon;
        DD2 = DD.*sign(Q);
        ID = abs(Q)>epsilon;
        E(ID) = DD2(ID);
    case 'Fro'
        E=(R2+mu*(X-X*S))/(2*beta+mu);
    case 'L21'
        alpha=beta/mu;
        G=X-X*S+R2/mu;
        G1 = sqrt(sum(G.^2,1));
        G1(G1==0) = alpha;
        G2 = (G1-alpha)./G1;
        E = G*diag((G1>alpha).*G2);        
end
error=E;