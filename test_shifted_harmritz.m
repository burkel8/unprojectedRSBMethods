which_matrix = 'largeLQCD';
N=40;
[A,n] = return_matrix(which_matrix,N);

s = 20;
base_shift = 0;
shift_increment = 100;
tol = 10e-02;
m = 20;
k= 10;

shift_recycle_method = 0;


fontsize = 13;
linewidth = 1;

X = rand(n,s);
B = rand(n,s);

shifts = zeros(1,s);
shift(1)= base_shift;

for i=2:s
shifts(i) = shifts(i-1) + shift_increment;
end

 Aclose = A + 0.1*sprand(A);
 K = rand(n,s);

 [Wtmp,Htmp] = block_Arnoldi(Aclose,K,m,s,n);
 
[P,~] = eigs(Htmp(1:m*s,1:m*s),k,'smallestabs');
U = Wtmp(:,1:m*s)*P;
C = A*U;
[C,R] = qr(C,0);
U = U/R;

[resid1] = unproj_rsbFOM(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method);
shift_recycle_method = 1;
[resid2] = unproj_rsbFOM(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method);
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



