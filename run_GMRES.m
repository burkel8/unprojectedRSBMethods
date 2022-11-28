%This run script solves a sequence of shifted linear systems using both
%sbFOM and unprojected rsbFOM. It counts the total number of vectors the
%matrix A is applied to to solve all systems in the sequence. For various
%values of the recycle subspace dimension, the results are plotted.

%User defined parameters
%  matrix - A string containing the name of the matrix to be used.
%           Current support for "smallLQCD", "sherman", "Poisson" and "chemical_potential"
matrix = 'smallLQCD';

%    N    - An integer used to specify the size of the Poisson matrix and
%    chemical potential matrix (value does not matter for other matrices)
N=80;

% num_systems - Number of systems to be solved in the sequence 
num_systems = 2;

% epsilon - parameter describing the strength of the pertubation between
% matrices in the sequence of problems. Default set to zero (special case
% where matrix does not change).
epsilon = 0;

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
k = [0,5,10];

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
[Aorig,n] = return_matrix(matrix,N,shift);

% Construct a num_shifts x s size matrix whose row i stores the 
% set of shifts of system i in the sequence.
shifts = zeros(num_systems,s);
shift = base_shift;

%populate entries of shifts matrix with appropriate shifts.
for i = 1:num_systems
  for j = 1:s
     shifts(i,j) = shift;
     shift = shift + shift_increment;
  end
end

numk = size(k,2); %numk stores the number of values of k we will test for.

%Initial solution approximation of zeros
X = zeros(n,s);

%matrix B storing right hand sides for each shift in the first system
B = rand(n,s);

%variables store final results
sbGMRES_results = zeros(1,numk);
unprojected_rsbGMRES_results = zeros(1,numk);
rsbGMRES_results = zeros(1,numk);

%We run the experiment for all values of k
for j = 1:numk
    fprintf("k = %d\n", k(j));

    %For each new experiment, we start off with 
    %  - the original matrix A
    %  - the total # of mat vecs set to 0
    %  - An empty recycling subspace U = []
    A = Aorig;
    total_sbGMRES_matvecs = 0;
    total_unprojected_rsbGMRES_matvecs = 0;
    total_rsbGMRES_matvecs = 0;
    U = [];

    %The U will need to be orthonormal for unprojected rsbGMRES, but not
    %for rsbGMRES. We do not want the U returned by unprojected rsbGMRES to
    %be used as the U for rsbGMRES. Thus denote them by two differnt names.
   
    U1= U; %recycling subspace for unprojected rsbGMRES
    U2 = U; %recycling subspace for rsbGMRES

%For each value of k, we solve the sequence of systems.
for i = 1:num_systems

fprintf("Solving system number %d\n", i);

%Solve each system using sbFOM and unprojected rsbFOM
[resid_sbGMRES, sbGMRES_matvecs] = sbGMRES(A,B,X,shifts(i,:),m,s,n,tol);
[resid_unprojected_rsbGMRES,U1,unprojected_rsbGMRES_matvecs] = unproj_rsbGMRES(A,B,X,shifts(i,:),m,k(j),s,n,tol,U1,shift_recycle_method);
[resid_rsbGMRES,U2,rsbGMRES_matvecs] = rsbGMRES(A,B,X,shifts(i,:),m,k(j),s,n,tol,U2,shift_recycle_method);

%Accumalate total number of vectors A is applied to throughout each
%algorithm.
total_sbGMRES_matvecs = total_sbGMRES_matvecs + sbGMRES_matvecs;
total_unprojected_rsbGMRES_matvecs = total_unprojected_rsbGMRES_matvecs + unprojected_rsbGMRES_matvecs;
total_rsbGMRES_matvecs = total_rsbGMRES_matvecs + rsbGMRES_matvecs;

%Generate new random right hand side for next system.
B = rand(n,s);

%Move to a new matrix
A = A + epsilon*sprand(A);
end

%Store total number of mat-vecs for both methods for each value of k
sbGMRES_results(j) = total_sbGMRES_matvecs;
unprojected_rsbGMRES_results(j) = total_unprojected_rsbGMRES_matvecs;
rsbGMRES_results(j) = total_rsbGMRES_matvecs;

end

%plot results
plot(k,sbGMRES_results,'LineWidth',linewidth);
hold on;
plot(k,unprojected_rsbGMRES_results,'LineWidth',linewidth);
hold on;
plot(k,rsbGMRES_results,'LineWidth',linewidth);
hold off;

title('sbGMRES vs. unproj rsbGMRES MAT-VECs','interpreter','latex','FontSize',fontsize)
xlabel('dim($\mathcal{U})$','interpreter','latex','FontSize',fontsize);
ylabel('MAT-VECs','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('sbGMRES','unprojected rsbGMRES','rsbGMRES','interpreter','latex');
set(lgd,'FontSize',fontsize);





