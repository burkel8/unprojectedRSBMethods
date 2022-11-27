%Function to solve eigenproblem required to compute k harmonic ritz vectors for the matrix A
% using the augmented Krylov subspace km(A,b) + U . 

%Input: m - dimension of Krylov subspace km(A,b)
%       k - number of eigenvectors to compute
%       G - augmented Hessenberg matrix G = [D 0; 0 Hbar];
%       What = [C Vm+1]
%       Vhat = [U Vm+1]

%Output: P - An m+k x k dimensional matrix storing the eigenvectors of the
%small eigenproblem in proposition 8.1 of preprint.
function [P] = harm_ritz_aug_krylov(m, k, G, What, Vhat)
    
    P = zeros(m+k,k);

    B = G'*(What'*What)*G;
    A = G'*(What'*Vhat);

    [harmVecs, harmVals] = eig(A,B);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),'descend');

    for i=1:k
        P(:,i) = harmVecs(:,iperm(i));
    end
end

