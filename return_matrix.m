%Function which returns the matrix to be used in each experiment.
%User can add a matrix here if they wish.
function [A,n] = return_matrix(which_matrix,N)

if which_matrix == "smallLQCD"
      load('smallLQCD_A1.mat');
      A=A1;
      [n,~] = size(A); 
      A = A - 0.65*speye(n);
      
   
elseif which_matrix == "largeLQCD"
      load('LQCD_8to4_0.mat');
     [n,~]=size(A);
     for i=1:n
       A(i,i) = A(i,i) - 0.59;
     end

elseif which_matrix == "poisson"
 
A = (N+1)^2*gallery('poisson',N);
n = N*N;

elseif which_matrix == "chemical_potential"
D2 = (N+1)^2*gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = D2 + 0*D1;
S = A;
[n,~]=size(S);
A = A - 19*speye(n);
else    
    fprintf('ERR: no matrix chosen!');
end

end


