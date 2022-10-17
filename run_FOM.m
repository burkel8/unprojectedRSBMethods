which_matrix = 'smallLQCD';
N=80;
[A,n] = return_matrix(which_matrix,N);


s = 5;
base_shift = 0;
shift_increment = 1;
tol = 10e-08;
m = 10;
k= [2,5,10];

shift_recycle_method = 0;


fontsize = 13;
linewidth = 1;

X = zeros(n,s);
B = rand(n,s);

shifts = zeros(1,s);
shift(1)= base_shift;

for i=2:s
shifts(i) = shifts(i-1) + shift_increment;
end

 Aclose = A + 0.1*sprand(A);
 K = rand(n,s);

 [Wtmp,Htmp] = block_Arnoldi(Aclose,K,m,s,n);
 [resid_sbFOM] = sbFOM(A,B,X,shifts,m,s,n,tol);

semilogy(resid_sbFOM,'-o','LineWidth',1);
hold on;

for i = 1:size(k,2)
     [P,~] = eigs(Htmp(1:m*s,1:m*s),k(i),'smallestabs');
U = Wtmp(:,1:m*s)*P;
C = A*U;
%[C,R] = qr(C,0);
%U = U/R;
[resid] = unproj_rsbFOM(A,B,X,shifts,m,k(i),s,n,tol,U,C,shift_recycle_method);
semilogy(resid,'-o','LineWidth',1);
hold on;

end
hold off;


hold off;
title('sbFOM vs. unprojected rsbFOM','interpreter','latex','FontSize',fontsize)
xlabel('Cycle Number','interpreter','latex','FontSize',fontsize);
ylabel('$\| \textbf{R} \|_{2}$','interpreter','latex','FontSize',fontsize);
grid on;
lgd = legend('sbFOM','un-proj rsbFOM $k=2$', 'un-proj rsbFOM $k=5$','un-proj rsbFOM $k=10$','interpreter','latex');
set(lgd,'FontSize',fontsize);



