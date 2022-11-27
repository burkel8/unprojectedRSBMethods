%% Note this code is used to examine how the spectrum of each matrix in the sequence
%% of problems changes for every value of epsilon.


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
num_quad_points = 1;   %number of quadrature points (add as many differnt points to this list)


matrix_eps = [0.0, 0.0001, 0.001, 0.01];  %parameter to determine how much the matrix changes.
num_systems = 30;

%Paramters for fontsize and line width in plots
fontsize = 13;
linewidth = 1;

%%%%%%%%%%%%%%    END USER INPUT HERE  %%%%%%%%%%%%%%%%%%%

[num_eps] = size(matrix_eps,2);
eigs_monitor = zeros(num_eps,num_systems);

%For each value of epsilon look at how the spectrum shifts for each
%new matrix
for jx = 1:num_eps

    [A,n] = return_matrix(which_matrix,N);

    for ix=1:num_systems

     A = A + matrix_eps(jx)*sprand(A);
     eigs_monitor(jx,ix) = real(eigs(A,1,'smallestreal'));
    
    end
end

%plot results
plot(eigs_monitor(1,:) ,'-s', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(eigs_monitor(2,:),'-v', 'LineWidth',1);
hold on;
plot(eigs_monitor(3,:),'-s', 'LineWidth',1,'MarkerSize', 8);
hold on;
plot(eigs_monitor(4,:),'-s', 'LineWidth',1);
hold off;

title('Minimum eigenvalue for each $\textbf{A}^{(i)}$','interpreter','latex', 'FontSize', fontsize)
xlabel('problem index','interpreter','latex', 'FontSize', fontsize);
ylabel('$ \mbox{min}(\lambda_{j}) $','interpreter','latex','FontSize',fontsize); 
grid on;
lgd = legend('$\epsilon = 0$', '$\epsilon = 0.0001$', '$\epsilon = 0.001$', '$\epsilon = 0.01$','interpreter','latex');
set(lgd,'FontSize',fontsize);

