%Function to apply the block Arnoldi process to build a block orthonormal 
% basis of the block Krylov subspace Km(A,B) 

%Input:
%       s - Block size
%       A - n x n matrix A
%       V - n x s normalized block vector such that R = VS for s x s matrix s
%       m - block Arnoldi cycle length
%       n - number of rows/cols of A

%Output
%       W - An n x (m+1)*s matrix storing a block orthonormal basis of Km(A,B)
%       H - Block upper Hessenberg matrix from block Arnoldi
%       nvm - number of vectors to which A was applied throughout the process.
function [W,H, nmv] = block_Arnoldi(A,V,m,s,n)

    nmv = 0; %accumalate number of vectors to which A is applied in this variable

    %preallocate W and H
    W = zeros(n,(m+1)*s);
    H = zeros((m+1)*s,m*s);

    %Assign V as first block element in W
    W(:,1:s) = V;
 
    %Modified Block Gram Schmidt block Arnoldi 
     for k = 1:s:m*s
        V = A*W(:,k:k+s-1);
        nmv = nmv + s;
        W(:,k+s:k+2*s-1) = V;
       
        for j=1:s:k
            H(j:j+s-1,k:k+s-1) = W(:,j:j+s-1)'*W(:,k+s:k+2*s-1);
            W(:,k+s:k+2*s-1) = W(:,k+s:k+2*s-1) - W(:,j:j+s-1)*H(j:j+s-1,k:k+s-1);
        end
        [W(:,k+s:k+2*s-1),H(k+s:k+2*s-1,k:k+s-1)] = qr(W(:,k+s:k+2*s-1),0);
     end
end