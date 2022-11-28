%Function which computes a recycling subspace U by solving the shift
%dependent harmonic Ritz problem outlined in proposition 7.1 of the
%preprint 

%Burke, L - Unprojected recycled block Krylov subspace methods for shifted systems  https://arxiv.org/abs/2209.06922

% Inputs:
%       H - block Hessenberg matrix
%       W - block Arnoldi basis
%       m - block Arnoldi cycle length
%       s - number of shifts
%       k - dimension of recycling subspace
%       U - recycling subspace
%       C - matrix such that C = A*U
%      shift - the shift sigma we wish to construct U for

% Output
%       U - recycling subspace.
function [U] = block_harm_ritz(H,W,m,s,k,U,C,shift)

eyebar = zeros((m+1)*s,m*s);
for i=1:m*s
eyebar(i,i) = 1;
end

%Construct eigenproblem
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

%Solve eigenproblem 
[harmVecs, harmVals] = eig(A1,A2);
harmVals = diag(harmVals);

%sort eigenvalues
[~,iperm] = sort(abs(harmVals),'descend');
for i=1:k
  P(:,i) = harmVecs(:,iperm(i));
end
U = What*G*P;
end


