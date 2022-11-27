%Function to apply Arnoldi algorithm to build basis for Krylov subspace Km(A,b)

%Input: A,b Matrix and vector used to build Km(A,b)
%     : n - dimension of A
%     : m - Arnolci cycle length
%     : re_orthog - argument to specify if the Arnoldi vectors should be
%     reorthogonalized. Set to 0 (no-reorthogonalization) or set to 1
%     (reorthogonalization).

%Output:  H,V - Hessenberg matrix and Arnoldi basis built from Arnoldi process. 

function [H,V] = arnoldi( A, b , n, m, re_orthog)

V = zeros(n,m+1);
H = zeros(m+1,m);

V(:,1) = b/norm(b);
for j=1:m

   V(:,j+1) = A*V(:,j);
  
   %modified GS
   for i=1:j
       H(i,j)=V(:,i)'*V(:,j+1);
       V(:,j+1)=V(:,j+1) - V(:,i)*H(i,j);
   end
    
    H(j+1,j) = norm(V(:,j+1));
    
    %optional reorthogonalization for stability 
    if re_orthog == 1
       for i=1:j
         mu = V(:,i)'*V(:,j+1);
         V(:,j+1) = V(:,j+1) - V(:,i)*mu;
         H(i,j) = H(i,j) + mu;
       end
        H(j+1,j) = norm(V(:,j+1));
    end
     V(:,j+1) = V(:,j+1)/H(j+1,j);
end
end
