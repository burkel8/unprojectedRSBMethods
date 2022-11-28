% Function which runs the unprojected recycled shifted block GMRES (unprojected rsbGMRES).

%Inputs 
%     A  - n x n Coefficient matrix from shifted linear system (A + sig*I)x(sig) = b(sig)
%     B  - n x s matrix whose columns are each right hand side b(sig)
%     X  - n x s matrix whose columns store an initial approximation for each sig
%    shifts - vector storing all shifts for the system.
%     m  - block Arnoldi cycle length
%     k  - recycling subspace dimension
%     s - number of shifts
%     n  - number of rows of A
%     tol - method convergence tolerance
%     U   - n x k matrix whos columns form a basis for recycling subspace
%    shift_recycle_method - parameter which (when set to 1) chooses to
%    build a recycling subspace appropriate to a new shift at each cycle.


%Output
%    resid - vector storing norm of block residual after each cycle.
%      U   - updated recycled subspace to be used in solving the next
%      system.
%      nvm - Total number of vectors the matrix A was applied to throughout

function [resid,U,nmv] = unproj_rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,shift_recycle_method)

nmv = 0; %variable to store total number of vectors matrix is applied to.
Y = zeros(m*s,s);

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

iter = 1;             %variable to keep track of iteration number
resid(iter) = norm(R);%assign norm of initial residual to first entry in output vector  
shift_monitor = 1;    %variable to store the index of the shift we wish to construct the recycling subspace
                      %for at each iteration. We will always construct U
                      %for the base shift unless the input variable
                      %shift_recycle_method is set to 1, in which case
                      %shift_monitor increments for each cycle.

%If we are on the first cycle of the first system with no previous
%recycling subspace to use, then run a cycle of sbFOM and use this cycle
% to construct a recycling subspace for the next cycle.
if(isempty(U))
   [r,X,R,H,W, sbGMRES_nmv] = sbGMRES_cycle(A,X,R,shifts,m,s,n);

   %Increment nmv by number of vectors A is applied to in sbFOM cycle.
   nmv = nmv + sbGMRES_nmv;
   iter = iter + 1;          

   %Add norm of updated residual in next element of residual vector.
   resid(iter) = r;

   %compute first augmentation subspace U using the eigenvalues of H.
   [P,~] = eigs(H(1:m*s,1:m*s),k,'smallestabs');
   [U] = W(:,1:m*s)*P;

   %Construct C (which costs an addition k MAT-Vecs)
   C = A*U;
   nmv = nmv + k;
else
   %If we already have a U from a previous system, then 
   % construct C for this new system (which costs an addition k MAT-Vecs).
   C = A*U;
   nmv = nmv + k;
end


while(norm(R) > tol)

%Compute normalized V to be used as starting block vector for block Arnoldi
% Store corresponding "R" factor, in s x s matrix S  
[V,~] = qr(R,0);

%Apply block Arnoldi process to build a basis for the block Krylov subspace Km(A,V)
[W,H,barnoldi_nmv] = block_Arnoldi(A,V,m,s,n);
nmv = nmv + barnoldi_nmv;

%precompute quantities.
Wp1TC = W'*C;
Wp1TU = W'*U;

Hsig = @(zx) (zx*eyebar + H);
CTU = @(zx) (C+zx*U)'*(C+zx*U);

%For each shift, solve problem on Krylov subspace and update solution and residual.
for i=1:s
sig = shifts(i);
Y(:,i) = (Hsig(sig)'*Hsig(sig) - Hsig(sig)'*(Wp1TC + sig*Wp1TU)*((CTU(sig))\((C+sig*U)'*W*Hsig(sig))))\...
    (Hsig(sig)'*W'*R(:,i) - Hsig(sig)'*(Wp1TC + sig*Wp1TU)*((CTU(sig))\((C+sig*U)'*R(:,i))));

X(:,i) = X(:,i) + W(:,1:m*s)*Y(:,i) + U*((CTU(sig)\((C+sig*U)'*(R(:,i) - W*Hsig(sig)*Y(:,i)))));

R(:,i) = R(:,i) - W*Hsig(sig)*Y(:,i) - (C+sig*U)*(((CTU(sig))\ (C+sig*U)'*(R(:,i) - W*Hsig(sig)*Y(:,i))));
end


%store updated residual norm in residual vector.
iter = iter + 1;
resid(iter) = norm(R);
shift = shifts(shift_monitor);

%Construct updated U by solving a harmonic ritz problem
[U] = block_harm_ritz(H,W,m,s,k,U,C,shift);

%Construct new C
C = A*U;
nmv = nmv + k; 

%Optionally select a recycling subspace suitable for a differnt shift on
%the next iteration.
if shift_recycle_method == 1
   shift_monitor = shift_monitor + 1;
   if shift_monitor == size(shifts,2)
       shift_monitor = 1;
   end
end

end