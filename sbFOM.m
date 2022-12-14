% Function which runs the shifted block Full Orthogonalization Method (sbFOM).

%Inputs 
%     A  - n x n Coefficient matrix from shifted linear system (A + sig*I)x(sig) = b(sig)
%     B  - n x s matrix whose columns are each right hand side b(sig)
%     X  - n x s matrix whose columns store an initial approximation for each sig
%    shifts - vector storing all shifts for the system.
%     m  - block Arnoldi cycle length
%     s  - number of shifts
%     n  - number of rows of A
%     tol - method convergence tolerance

%Output
%    resid - vector storing norm of block residual after each cycle.
%      nvm - Total number of vectors the matrix A was applied to throughout

function [resid, nmv] = sbFOM(A,B,X,shifts,m,s,n,tol)

nmv = 0; %variable to accumalate total number of vectors A is applied to throughout sbFOM
Y = zeros(m*s,s);
E = zeros(m*s,s);

%create (j+1)s x js matrix of zeros with identity matrix on first js x js rows
eyebar = zeros(m*s+s,m*s);
for i=1:m*s
eyebar(i,i)=1;
end

%allocate memory for block residual vector
R=zeros(n,s);

%Construct initial residual for each shift in the system.
for i=1:s
R(:,i) = B(:,i) - (A + shifts(i)*speye(n))*X(:,i);
nmv = nmv + 1;
end

iter = 1;  %variable to keep track of iteration number

%assign norm of initial residual to first entry in output vector  
resid(iter) = norm(R);

%Perform enough cycles of the method until the residual norm is descresed below the tolerance level.                      
while(norm(R) > tol)

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

Y(:,i) = (H(1:m*s,1:m*s) + sig*speye(m*s))\E(:,i);
X(:,i) = X(:,i) + W(:,1:m*s)*Y(:,i);
R(:,i) = R(:,i) - W*Hsig(sig)*Y(:,i);
end

%store updated residual norm in residual vector.
iter = iter + 1;
resid(iter) = norm(R);

end
end