%This program tests all three versions of rFOM2 as a recycle method by approximating
% a sequence of f(A)b problems and recycling between them.

%Choose parameters for program

%%First choose the matrix
% Possible options are
%-- A small lattice QCD matrix of size 3072x3072 ("smallLQCD")
%-- A poisson matrix of size N*N x N*N (user specifies N) ("poisson")
%-- A chemical potential matrix of size N*N x N*N (user specifies N) ("chemical_potantial")
matrix = "smallLQCD";   

%%Choose the function . Available options are
% -- inverse function ("inverse")
% -- invSqrt function ("invSqrt")
% -- log function ("log")
% -- square root function ("sqrt")
problem = 'inverse';

%Shift the matrix by some multiple of the identity matrix. Do this to
%ensure the spectrum of the matrix remains positive.
%To see the full benifits of recycling choose the shift such that the
%the smallest eigenvalue is close to the origin. Some suggestions for the 
%Hermitian and non-Hermitian QCD matrices are given below. Care should be
%taken when changing these values.
if strncmp(matrix,"smallLQCD",20) == 1
   shift =  0.065;
elseif strncmp(matrix,"hermitian_QCD",20) == 1
   shift = -7.7;
else % 0 to be used for all other matrices
   shift = 0;
end

%Suggested number of quadrature points for invSqrt and inverse
%functions are given below. 
if strncmp(problem,"inverse",20) == 1
  num_quad_points = 3000;   
elseif strncmp(problem,"invSqrt",20) == 1
  num_quad_points = 30; 
else
  num_quad_points = 30; 
end

%% Parameters of solve
m = 40;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 50;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)


matrix_eps = 0.0;  %parameter to determine how much the matrix changes.
num_systems = 7;


%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%
e1 = zeros(m,1);
e1(1)=1;

%vectors to store results
err_arnoldi = zeros(1,num_systems);
err_quad_arnoldi = zeros(1,num_systems);
err_rFOM_v1 = zeros(1,num_systems);
eigs_monitor = zeros(1,num_systems);

%store matrix and function in appropriate vectors
[f_scalar, f_matrix] = return_function(problem);
[A,n] = return_matrix(matrix,N,shift);

%create vector
b = rand(n,1);
b = b/norm(b);

%Compute exact solution of first problem.
x = f_matrix(A,b);

%Run Arnoldi on a close matrix to generate first augmentation subspace U
Aclose = A + 0.001*sprand(A);
g = rand(n,1);
[Hc,Vc] = arnoldi( Aclose, g , n,m, 1);

%Construct initial U
%[P] = harmonic_ritz(Hc,m,k);
[P,~] = eigs(Hc(1:m,1:m),k,'smallestabs');
U = Vc(:,1:m)*P;
C = A*U;

%Create sequence of problems and approximate f(A)b for each.
for ix=1:num_systems

    %Variable to monitor the smallest real eigenvalue of each matrix just
    %to ensure the spectrum does not shift ti the negative real part of the
    %complex plane. Line not essential.
    eigs_monitor(ix) = real(eigs(A,1,'smallestreal'));

    fprintf("\n Approximating f(A)b # %d ... \n\n",ix);

    %Run Arnoldi to build basis V (of Krylov subspace) and Hessengerg matrix H
    [H,V] = arnoldi( A, b , n,m, 1);
  
    %Standard Arnoldi Approximation
    arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
    err_arnoldi(ix) = norm(x - arnoldi_approx);

    %Quadrature Arnoldi approximation

    %For the inverse square root, use special quadrature, else use
    %trapezoidal rule
    if strncmp(problem,"invSqrt",20) == 1
     quad_arnoldi_Approx = quad_arnoldi_invSqrt(V,H,m,num_quad_points);
    else 
      quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad_points,f_scalar);
    end
    err_quad_arnoldi(ix) = norm(x - quad_arnoldi_Approx);

    %r(FOM)2 v1
    %For the inverse square root, use special quadrature, else use
    %trapezoidal rule
    if strncmp(problem,"invSqrt",20) == 1
    [rFOM_v1_approx] = rFOM2_v1_invSqrt(b,V,H,m,k,U,C,num_quad_points);
    else 
     [rFOM_v1_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad_points,f_scalar);
    end
    err_rFOM_v1(ix) = norm(x - rFOM_v1_approx);

    fprintf("\n... DONE\n");

     %Construct G
         [U,D] = scale_cols_of_U(U,k);
         Vhat = [U V(:,1:m)];
         What = [C V(:,1:m+1)];
         G = zeros(m+1+k,m+k);
         G(1:k,1:k) = D;
         G(k+1:m+1+k,k+1:m+k) = H;
    
    %Create new problem in sequence and compute its exact solution
    b = rand(n,1);
    b = b/norm(b);

    %If we are working with a Hermitian QCD matrix make a random
    %pertubation to ensure result is still Hermitian.
    if strncmp(matrix,"hermetian_QCD",20) == 1
    randVec = eps*rand(n,1);
    A = A + toeplitz(randVec);
    else
    %if matrix is non-Hermitian, no need to make a symmetric pertubation.
    A = A + matrix_eps*sprand(A);
    end
   
    %Compute exact solution of new problem.
    x = f_matrix(A,b);

    %Compute new U for next problem in the sequence
    [P] = harm_ritz_aug_krylov(m,k,G,What,Vhat);
    U = Vhat*P;
    C = A*U;
end

%Plot Results
semilogy(err_arnoldi/norm(x) ,'-s', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
semilogy(err_quad_arnoldi/norm(x) ,'-o', 'LineWidth',1);
hold on;
semilogy(err_rFOM_v1/norm(x),'-v', 'LineWidth',1);
hold off;

title('sign($\textbf{A}$)\textbf{b} - error vs. problem index','interpreter','latex', 'FontSize', fontsize)
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('$\| f(\textbf{A})\textbf{b} - \textbf{x}_{m} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('Arnoldi','Arnoldi (q)','rFOM$^{2}$','interpreter','latex');
set(lgd,'FontSize',fontsize);


