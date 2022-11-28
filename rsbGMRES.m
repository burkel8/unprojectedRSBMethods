%An implementation of the rsbGMRES algorithm from 

%Soodhalter, Kirk M. "Block Krylov subspace recycling for shifted systems with 
% unrelated right-hand sides."SIAM Journal on Scientific Computing 38.1 (2016): A302-A324.

%This code has been adapted from the code associated to the above
%paper. (Available at https://zenodo.org/record/56157#.Y4TTJnbP1PY)

%Input: 
% Function which runs the unprojected recycled shifted block Full
% Orthogonalization Method (unprojected rsbFOM).

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

function [resid,U,nmv] = rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,shift_recycle_method)

nmv = 0;

%create (j+1)s x js matrix of zeros with identity matrix on first js x js rows
eyebar = zeros(m*s+s,m*s);
for i=1:m*s
eyebar(i,i)=1;
end

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
   
   C = A*U;
   [C,RR] = qr(C,0);
   U = U/RR;
   nmv = nmv + k;
else
   %If we already have a U from a previous system, then 
   % construct C for this new system (which costs an addition k MAT-Vecs).
   C = A*U;
   nmv = nmv + k;
end

%Perform enough iterations of the method until the residual norm is small enough
while(norm(R) > tol)


 %perform initial residual (and solution) projections to ensure residuals
 %are orthogonal to C = A*U
 CTU = C'*U;
 CTR = C'*R;
 for i = 1:length(shifts)
     rProj = ( eye(size(U,2)) + shifts(i)*CTU ) \ CTR(:,i);
     X(:,i) = X(:,i) + U*rProj;
     R(:,i) = R(:,i) - C*rProj - U*(shifts(i)*rProj);
 end

%Built basis for block Krylov subspace Km(A,R)
[V,S0] = qr(R,0);

%We use a differnt block Arnoldi function to the one we used with the unprojected
%methods. This rsbGMRES Arnoldi performs an additional orthogonalization of
%the Arnoldi block vectors with respect to C.
[W,H,B,nmv_barnoldi] = rsb_GMRES_block_arnoldi(A,V,m,C);

nmv = nmv + nmv_barnoldi;

 WTU = W(:,1:(m+1)*s)'*U;
 Uhat = U - C*CTU - W(:,1:(m+1)*s)*WTU;
 Vhat = orth(Uhat);
 VTU = Vhat'*U;
 kr = size(Vhat,2);
 RES = zeros((m+1)*s+k+kr,s);
 RHS = eye((m+1)*s+k+kr,s);
 RHS(1:s,:) = S0;
 RHS(end-kr+1:end,:) = Vhat'*R;
 Y = zeros(m*s+k,s);
        
 HH = @(sigma) H(1:(m+1)*s,1:m*s) + sigma*eye((m+1)*s,m*s);
 G = @(sigma) [HH(sigma) sigma*WTU; B(:,1:m*s) eye(k)+sigma*CTU; zeros(kr,m*s) sigma*VTU];
        
        for p=1:s
            if shifts(p) ~= 0
                Gsig = G(shifts(p));
                Y(:,p) = Gsig\RHS(:,p);
                RES(:,p) = RHS(:,p) - Gsig*Y(:,p);
            else
                rhstemp = zeros((m+1)*s,1);
                rhstemp(1:s) = S0(:,p);
                ytemp = H(1:(m+1)*s,1:m*s)\rhstemp;
                Y(1:m*s,p) = ytemp;
                Y(m*s+1:end,p) = - B(:,1:m*s)*ytemp;
                restemp = rhstemp - H(1:(m+1)*s,1:m*s)*ytemp;
                RES(1:(m+1)*s,p) = restemp;
            end
        end

    % Compute approximation at end of cycle
    X = X + W(:,1:m*s)*Y(1:m*s,:);
    X = X + U*Y(m*s+1:end,:);
    R = W(:,1:(m+1)*s)*RES(1:(m+1)*s,:);
    R = R + C*RES((m+1)*s+1:(m+1)*s+k,:);
    R = R + Vhat*RES((m+1)*s+k+1:end,:);

iter = iter + 1;
resid(iter) = norm(R);

%compute new recycle subspace for next cycle.
shift = shifts(shift_monitor);
[U] = block_harm_ritz(H,W,m,s,k,U,C,shift);
C = A*U;
[C,RR] = qr(C,0);
U = U/RR;
nmv = nmv + k;

 if shift_recycle_method == 1
    shift_monitor = shift_monitor + 1;
    if shift_monitor == size(shifts,2) 
        shift_monitor = 1;
    end
 end
end
end





