%An implementation of the block Arnoldi process for the rsbGMRES algorithm from

%Soodhalter, Kirk M. "Block Krylov subspace recycling for shifted systems with 
% unrelated right-hand sides."SIAM Journal on Scientific Computing 38.1 (2016): A302-A324.

%In this Arnoldi process, the block Arnoldi vectors are also orthogonalized
%against C.

%This code has been adapted from the code associated to the above
%paper. (Available at https://zenodo.org/record/56157#.Y4TTJnbP1PY)

%Input:
%       A - n x n matrix A
%       V - n x s normalized block vector such that R = VS for s x s matrix s
%       m - block Arnoldi cycle length
%       C - matrix C such that C = A*U

%Output
%       W - An n x (m+1)*s matrix storing a block orthonormal basis of Km(A,B)
%       H - Block upper Hessenberg matrix from block Arnoldi
%       B - k x m block upper Hessenberg matrix storing orthogonalization
%       coefficients against C
%       nvm - number of vectors to which A was applied throughout the process.
function [W,H,B,nmv] = rsb_GMRES_block_arnoldi(A,V,m,C)

    nmv = 0;  %variable to store number of vectors matrix is applied to.
    [n,s] = size(V);
    k = size(C,2);
    H = zeros((m+1)*s,m*s);
    W = zeros(n,(m+1)*s);
    B = zeros(k,m*s);
    
    W(:,1:s) = V;

    iter = 0;
    %%%%%%%%%%%%%%%%%%%%
    % Iterate by block %
    %%%%%%%%%%%%%%%%%%%%
    for pp = 1:s:m*s
        iter = iter + 1;
        
        V = A*W(:,pp:pp+s-1);
        nmv = nmv + s;
        W(:,pp+s:pp+2*s-1) = V;
        
        %%%%%%%%%%%%%%%%%%%%%
        % Orthog against C  %
        % Store coeffs in B %
        %%%%%%%%%%%%%%%%%%%%%
        B(:,pp:pp+s-1) = C(:,1:k)'*W(:,pp+s:pp+2*s-1);
        W(:,pp+s:pp+2*s-1) = W(:,pp+s:pp+2*s-1) - C(:,1:k)*B(:,pp:pp+s-1);
        
        %%%%%%%%%%%%%%%%%
        % Block Arnoldi %
        % Process       %
        %%%%%%%%%%%%%%%%%
        for j=1:s:pp
            H(j:j+s-1,pp:pp+s-1) = W(:,j:j+s-1)'*W(:,pp+s:pp+2*s-1);
            W(:,pp+s:pp+2*s-1) = W(:,pp+s:pp+2*s-1) - W(:,j:j+s-1)*H(j:j+s-1,pp:pp+s-1);
        end
        [W(:,pp+s:pp+2*s-1),H(pp+s:pp+2*s-1,pp:pp+s-1)] = qr(W(:,pp+s:pp+2*s-1),0);
    end
end
