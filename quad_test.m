%This program tests all three versions of r(FOM)2 as an augmented Krylov 
%subspace method. User is required to choose parameters for program below.

%%First choose the matrix 
% possible options are
%-- A small lattice QCD matrix of size 3072x3072 ("smallLQCD")
%-- A poisson matrix of size N*N x N*N (user specifies N) ("poisson")
%-- A chemical potential matrix of size N*N x N*N (user specifies N) ("chemical_potantial")
matrix = "hermitian_QCD";   

%%Choose the function 
% Possible options are
% -- inverse function ("inverse")
% -- invSqrt function ("invSqrt")
% -- log function ("log")
% -- square root function ("sqrt")
problem = 'invSqrt';

m = 40;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 100;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)

         
%Shift the matrix by some multiple of the identity matrix. Do this to
%ensure the spectrum of the matrix remains positive.
%To see the full benifits of recycling choose the shift such that
%the smallest eigenvalue is close to the origin. Some suggestions for the 
%Hermitian and non-Hermitian QCD matrices are given below. Care should be
%taken when changing these values.
if strncmp(matrix,"smallLQCD",20) == 1
   shift =  0.65;
elseif strncmp(matrix,"hermitian_QCD",20) == 1
   shift = -7.7;
else % 0 to be used for all other matrices
   shift = 0;
end

%Suggested number of quadrature points for invSqrt and inverse
%functions are given below. 
if strncmp(problem,"inverse",20) == 1
   num_quad = [100,200,300,500,900];
elseif strncmp(problem,"invSqrt",20) == 1
 num_quad = [1,3,5,7,8,10,15,20,30,40,50]; 
 else
  num_quad = [10000,20000]; 
end

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;
%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%
%Store matrix and function in appropriate variables.
[A,n] = return_matrix(matrix,N,shift);
[f_scalar, f_matrix] = return_function(problem);

% Define vector
b = rand(n,1);
b = b/norm(b);

num_tests = size(num_quad,2); 

%compute exact solution
exact = f_matrix(A,b);
e1 = zeros(m,1); e1(1)=1;

%Create a augmentation subspace directly from A using eigs
[U,~] = eigs(A,k,'smallestabs');
C = A*U;

%vectors to store results of each approximation.
err_arnoldi = zeros(1,num_tests);
err_quad_arnoldi = zeros(1,num_tests);
err_rFOM_v1 = zeros(1,num_tests);
err_rFOM_v2 = zeros(1,num_tests);
err_rFOM_v3 = zeros(1,num_tests);

%Run Arnoldi
[H,V] = arnoldi( A, b , n,m, 1);

%Repeat each experiment for differnt numbers of quadrature points
for i=1:num_tests

% Compute approximation using standard Arnoldi approximation
arnoldi_approx = norm(b)*V(:,1:m)*f_matrix(H(1:m,1:m),e1);
err_arnoldi(i) = norm(exact - arnoldi_approx);

% Compute approximation using quadrature Arnoldi approximation
%For the inverse square root, use special quadrature, else use trapezoidal rule
 if strncmp(problem,"invSqrt",20) == 1
     quad_arnoldi_Approx = quad_arnoldi_invSqrt(V,H,m,num_quad(i));
    else 
      quad_arnoldi_Approx = quad_arnoldi(b,V,H,m,num_quad(i),f_scalar);
 end
 err_quad_arnoldi(i) = norm(exact - quad_arnoldi_Approx);

%Compute approximation r(FOM)2 version 1
%quadrature Arnoldi approximation
%For the inverse square root, use special quadrature, else use trapezoidal rule
 if strncmp(problem,"invSqrt",20) == 1
  [rFOM_v1_approx] = rFOM2_v1_invSqrt(b,V,H,m,k,U,C,num_quad(i));
 else 
  [rFOM_v1_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad(i),f_scalar);
 end
   err_rFOM_v1(i) = norm(exact - rFOM_v1_approx);

%Compute approximation using r(FOM)2 version 2.
%quadrature Arnoldi approximation
%For the inverse square root, use special quadrature, else use trapezoidal rule
if strncmp(problem,"invSqrt",20) == 1
    [rFOM_v2_approx] = rFOM2_v2_invSqrt(b,V,H,m,k,U,C,num_quad(i));
    else 
     [rFOM_v2_approx] = rFOM2_v2(b,V,H,m,k,U,C,num_quad(i), f_scalar);
end
  err_rFOM_v2(i) = norm(exact - rFOM_v2_approx);

%Compute approximation using r(FOM)2 version 3
%quadrature Arnoldi approximation
%For the inverse square root, use special quadrature, else use trapezoidal rule
 if strncmp(problem,"invSqrt",20) == 1
    [rFOM_v3_approx] = rFOM2_v3_invSqrt(b,V,H,m,k,U,C,num_quad(i),f_matrix);
    else 
     [rFOM_v3_approx] = rFOM2_v3(b,V,H,m,k,U,C,num_quad(i), f_scalar, f_matrix);
  end
err_rFOM_v3(i) = norm(exact - rFOM_v3_approx);

end

%plot results
x=1:1:num_tests;
points = num_quad(x);
semilogy(points,err_arnoldi,'-s','LineWidth',1);
hold on;
semilogy(points,err_quad_arnoldi,'-o','LineWidth',1);
hold on;
semilogy(points,err_rFOM_v1,'-v','LineWidth',1 ,'MarkerSize', 9);
hold on;
semilogy(points,err_rFOM_v2,'-o','LineWidth',1);
hold on;
semilogy(points,err_rFOM_v3,'-s','LineWidth',1);
hold off;
title(' Approximation Error ','interpreter','latex','FontSize',fontsize)
xlabel('number of quad nodes','interpreter','latex','FontSize',fontsize);
ylabel('$\| f(\textbf{A})\textbf{b} - \tilde{f}_{i} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('Arnoldi','Arnoldi (q)','rFOM$^{2}$ $\tilde{f}_{1}$','rFOM$^{2}$ $\tilde{f}_{2}$', 'rFOM$^{2}$ $\tilde{f}_{3}$','interpreter','latex');
set(lgd,'FontSize',fontsize);
xticks(points)

