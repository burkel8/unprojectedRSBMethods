%This run script compares the convergence for a single system of the unprojected rsbGMRES 
% algorithm to the rsbGMRES algorithm of [Soodhalter,2016,SISC]. 

%User defined parameters
%  matrix - A string containing the name of the matrix to be used.
%           Current support for "smallLQCD", "sherman", "Poisson" and "chemical_potential"
matrix = 'smallLQCD';

%    N    - An integer used to specify the size of the Poisson matrix and
%    chemical potential matrix (value does not matter for other matrices)
N=80;

% s - number of shifts for each system (We assume all systems in the sequence have 
% same number of shifts).
s = 1;

% base_shift - the value of the first shift in the first system
base_shift = 0;

% shift_increment - the set of shifts for each system are a distance of
% shift_increment apart.
shift_increment = 0.1;

% tol - the residual norm convergence tolerance for the method 
tol = 10e-5;

% m - The block Arnolci cycle length
m = 20;

% k - vector containing all values of recycle subspace dimension we wish to test
k = 10;

% shift_recycle_method - a variable which when set to 1 constructs a
% recycle subspace corresponding to a new shift at each iteration of the
% method. By default, we leave this variable turned off, meaning we only
% recycle for the shift sig = 0.
shift_recycle_method = 0;

%parameters for plots.
fontsize = 13;
linewidth = 1;

%We may need to shift the spectrum of the input matrices. Some default
%options are set below for differnt matrices.
if strcmp(matrix,"sherman")==1
   shift = -190;
elseif strcmp(matrix,"smallLQCD")==1
  shift = 0.65;
end

%construct matrix A
[A,n] = return_matrix(matrix,N,shift);

% Construct a num_shifts x s size matrix whose row i stores the 
% set of shifts of system i in the sequence.
shifts = zeros(1,s);
shift = base_shift;

%populate entries of shifts vector with appropriate shifts.
  for j = 1:s
     shifts(1,j) = shift;
     shift = shift + shift_increment;
  end

%Initial solution approximation of zeros
X = zeros(n,s);

%Initial set of random vectors for each shift
B = rand(n,s);

%Initial recycling subspace is empty.
U = [];

%Call rsbGMRES
[rsbGMRES_resid] = rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,shift_recycle_method);

%Call unprojected rsbGMRES
[unproj_rsbGMRES_resid] = unproj_rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,shift_recycle_method);

%plot convergence curves.
semilogy(rsbGMRES_resid,'-o','LineWidth',1);
hold on;
semilogy(unproj_rsbGMRES_resid,'-o','LineWidth',1);
hold off
title('rsbGMRES  vs. unprojected rsbGMRES','interpreter','latex','FontSize',fontsize)
xlabel('Cycle Number','interpreter','latex','FontSize',fontsize);
ylabel('$\| \textbf{R} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('rsbGMRES ','un-proj rsbGMRES ','interpreter','latex');
set(lgd,'FontSize',fontsize);



