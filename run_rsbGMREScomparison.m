which_matrix = 'smallLQCD';
N=80;
[A,n] = return_matrix(which_matrix,N);

s = 3;
base_shift = 0;
shift_increment = 0.01;
tol = 10e-10;
m = 30;
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

[rsbGMRES_resid] = rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method);
[unproj_rsbGMRES_resid] = unproj_rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method);
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



