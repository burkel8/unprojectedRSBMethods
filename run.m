load('smallLQCD_A1.mat');
A=A1;
[n,~] = size(A);
b = rand(n,1);

s = 10;
base_shift = 0;
shift_increment = 0.1;
tol = 1.0e-09;
m = 20;
k=10;

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
 [P,~] = eigs(Htmp(1:m*s,1:m*s),k,'smallestabs');
U = Wtmp(:,1:m*s)*P;


[U,~] = eigs(A,k,'smallestabs');
C = A*U;
[C,R] = qr(C,0);
U = U/R;

[resid] = unproj_rsbFOM(A,B,X,shifts,m,k,s,n,tol,U,C);
semilogy(resid)
grid on;

