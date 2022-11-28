% Function which runs a single cycle of the shifted block GMRES (sbFOM).

%Inputs 
%     A  - n x n Coefficient matrix from shifted linear system (A + sig*I)x(sig) = b(sig)
%     X  - n x s matrix whose columns store an initial approximation for each sig
%     R  - n x s matrix storing initial block residual
%    shifts - vector storing all shifts for the system.
%     m  - block Arnoldi cycle length
%     s  - number of shifts
%     n  - number of rows of A

%Output
%       r  - block residual norm after the cycle
%       X  -  updated solution after the cycle.
%       R  - updated residual after the cycle.
%       H  - ms x ms Block upper Hessenberg matrix from block Arnoldi process
%       W  - n x (m+1)s matrix storing block orthonormal basis of block
%       Krylov subspace K(A,R).
%      nvm - Total number of vectors the matrix A was applied to throughout


function [r,X,R,H,W,nmv] = sbGMRES_cycle(A,X,R,shifts,m,s,n)

nmv = 0;
Y = zeros(m*s,s);
E = zeros(m*s+s,s);

%create (j+1)s x js matrix of zeros with identity matrix on first js x js rows
eyebar = zeros(m*s+s,m*s);
for i=1:m*s
eyebar(i,i)=1;
end

%Compute normalized V to be used as starting block vector for block Arnoldi
% Store corresponding "R" factor, in s x s matrix S  
[V,S] = qr(R,0);

%Apply block Arnoldi process to build a basis for the block Krylov subspace Km(A,V)
[W,H, barnoldi_nmv] = block_Arnoldi(A,V,m,s,n);
nmv = nmv + barnoldi_nmv;
E(1:s,1:s) = S;

Hsig = @(zx) H + zx*eyebar;

%For each shift, solve small problem on Krylov subspace. Update solution
%and residual accordingly.
for i=1:s
sig = shifts(i);

Y(:,i) = (H + sig*eyebar)\E(:,i);
X(:,i) = X(:,i) + W(:,1:m*s)*Y(:,i);
R(:,i) = R(:,i) - W*Hsig(sig)*Y(:,i);
end

r = norm(R);