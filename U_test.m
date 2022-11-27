%% Note this code is used to test the accuracy of U as a recycle subspace
%% As we move through the sequence of f(A)b problems. The accuracy of U does not
%% depend on the actual evaluation of f(A)b using our methods, so we do not include
%% these computations in the code for speed.

%% Step 1: Choose parameters for program

%%First choose the matrix possible options are
%-- A small lattice QCD matrix of size 3072x3072 ("smallLQCD")
%-- A poisson matrix of size N*N x N*N (user specifies N) ("poisson")
%-- A chemical potential matrix of size N*N x N*N (user specifies N) ("chemical_potantial")
which_matrix = "smallLQCD";   

%%Choose the function . Available options are
% -- inverse function ("inverse")
% -- Sign function ("sign")
% -- log function ("log")
% -- square root function ("sqrt")
problem = 'inverse';


%% Parameters of solve
m = 40;  %Arnoldi cycle length
k = 20;  %recycle space dimension
N = 50;  %Parameter for Poisson and chemical potential matrix (value 
         %does not matter for other matrices)
num_quad_points = 1;   %This experiment does not depend on the number of quad points. Choose 1 for speed


matrix_eps = [0.0, 0.0001, 0.001, 0.01];  %parameter to determine how much the matrix changes.
num_systems = 30;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%

[num_eps] = size(matrix_eps,2);
err_recycle_space = zeros(num_eps,num_systems);

%For each value of epsilon, approximate a sequence of f(A)b applications
%and generate a new U for each problem using the previous problem. Store
%the error of U as an eigenvector approximation for epsilon index i and
%problem index j in the entry (i,j) of the matrix err_recycle_space.
for jx = 1:num_eps

[S1,n] = return_matrix(which_matrix,N);

b = rand(n,1);
b = b/norm(b);

%% Run Arnoldi on a close matrix
Sclose = S1 + 0.001*sprand(S1);
g = rand(n,1);
[Hc,Vc] = arnoldi( Sclose, g , n,m, 1);
[P] = harmonic_ritz(Hc,m,k);
U = Vc(:,1:m)*P;
C = S1*U;
S = S1;

for ix=1:num_systems
  
    fprintf("\nSOLVING FOR SYSTEM # %d ... \n\n",ix);
    [H,V] = arnoldi( S, b , n,m, 1);
    [U,D] = scale_cols_of_U(U,k);
    Vhat = [U V(:,1:m)];
    What = [C V(:,1:m+1)];
    G = zeros(m+1+k,m+k);
    G(1:k,1:k) = D;
    G(k+1:m+1+k,k+1:m+k) = H;
        
    %Change to new problem
    b = rand(n,1);
    b = b/norm(b);
    S = S + matrix_eps(jx)*sprand(S);
    [V,~] = eigs(S1,k,'smallestabs');


    %Compute new U using previous Arnoldi
    [P] = harm_ritz_aug_krylov(m,k,G,What,Vhat);
    U = Vhat*P;
    [U,R] = qr(U,0);
    C = S*U;
    %[C,R] = qr(S*U,0);
    %U = U/R;

    %measure error
    %theta = (U'*U)\(U'*S*U);
    %eigs_res_norm = norm(U*theta - S*U)/norm(U);
    %fprintf("avg_res_norm = %f\n", eigs_res_norm);

    eigs_res_norm = subspace(U,V);

    %store result
    err_recycle_space(jx,ix) = eigs_res_norm;
    
end
end

%Plot results
semilogy(err_recycle_space(1,:) ,'-s', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
semilogy(err_recycle_space(2,:),'-v', 'LineWidth',1);
hold on;
semilogy(err_recycle_space(3,:),'-s', 'LineWidth',1,'MarkerSize', 8);
hold on;
semilogy(err_recycle_space(4,:),'-s', 'LineWidth',1);
hold off;

title('Accuracy of $\mathcal{U}$ as an eigenvector approximation','interpreter','latex', 'FontSize', fontsize)
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('$\theta(\textbf{U},\textbf{Z})$','interpreter','latex','FontSize',fontsize); 
grid on;
lgd = legend('$\epsilon = 0$', '$\epsilon = 0.0001$', '$\epsilon = 0.001$', '$\epsilon = 0.01$','interpreter','latex');
set(lgd,'FontSize',fontsize);

