%This run script tests the effects on covergence of generating a recycling
%subspace U corresponding to a new shift at each cycle of the unprojected rsbFOM.

%User defined parameters
%  matrix - A string containing the name of the matrix to be used.
%           Current support for "smallLQCD", "sherman", "Poisson" and "chemical_potential"
matrix = 'smallLQCD';

%    N    - An integer used to specify the size of the Poisson matrix and
%    chemical potential matrix (value does not matter for other matrices)
N=80;

% s - number of shifts for each system (We assume all systems in the sequence have 
% same number of shifts).
s = 10;

% base_shift - the value of the first shift in the first system
base_shift = 0;

% shift_increment - the set of shifts for each system are a distance of
% shift_increment apart.
shift_increment = 100;

% tol - the residual norm convergence tolerance for the method 
tol = 10e-5;

% m - The block Arnolci cycle length
m = 20;

% k - vector containing all values of recycle subspace dimension we wish to test
k = 10;

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

%start with initial empty recycling subspace.
U = [];

% shift_recycle_method - a variable which when set to 1 constructs a
% recycle subspace corresponding to a new shift at each iteration of the
% method. We first run the function with this variable set to 0.
shift_recycle_method = 0;
[resid1] = unproj_rsbFOM(A,B,X,shifts,m,k,s,n,tol,U,shift_recycle_method);
%Then run with it set to 1.
shift_recycle_method = 1;
[resid2] = unproj_rsbFOM(A,B,X,shifts,m,k,s,n,tol,U,shift_recycle_method);

%plot results
semilogy(resid1,'-o','LineWidth',1);
hold on;
semilogy(resid2,'-o','LineWidth',1);
hold off
title('unproj-rsbGMRES','interpreter','latex','FontSize',fontsize)
xlabel('Cycle Number','interpreter','latex','FontSize',fontsize);
ylabel('$\| \textbf{R} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('unproj rsbGMRES m1','un-proj rsbGMRES m2','interpreter','latex');
set(lgd,'FontSize',fontsize);



