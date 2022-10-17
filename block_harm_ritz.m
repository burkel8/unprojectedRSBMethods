function [U] = block_harm_ritz(H,W,m,s,k,U,C,shift)


eyebar = zeros((m+1)*s,m*s);
for i=1:m*s
eyebar(i,i) = 1;
end
n = size(W,1); 
Rz = zeros(n,k+m*s);
P = zeros(k+m*s,k);
G = zeros(k+(m+1)*s,k+m*s);
G(1:k,1:k) = speye(k);
G(k+1:(m+1)*s+k,k+1:m*s+k) = H + shift*eyebar;
Rz(:,1:k) = shift*U;
Vhat = [U W(:,1:m*s)];
What = [C W];
A1 = (What*G + Rz)'*Vhat;
A2 = (What*G + Rz)'*(What*G+Rz);

 [harmVecs, harmVals] = eig(A1,A2);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),'descend');

    for i=1:k
        P(:,i) = harmVecs(:,iperm(i));
    end

    U = What*G*P;
end


