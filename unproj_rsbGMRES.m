function [resid] = unproj_rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method)

Y = zeros(m*s,s);
eyebar = zeros(m*s+s,m*s);
for i=1:m*s
eyebar(i,i)=1;
end

R=zeros(n,s);

for i=1:s
R(:,i) = B(:,i) - (A + shifts(i)*speye(n))*X(:,i);
end

iter = 1;
resid(iter) = norm(R);

shift_monitor = 1;


while(norm(R) > tol)

[W,H] = block_Arnoldi(A,R,m,s,n);


Wp1TC = W'*C;
Wp1TU = W'*U;

Hsig = @(zx) (zx*eyebar + H);
CTU = @(zx) (C+zx*U)'*(C+zx*U);

for i=1:s
sig = shifts(i);
Y(:,i) = (Hsig(sig)'*Hsig(sig) - Hsig(sig)'*(Wp1TC + sig*Wp1TU)*((CTU(sig))\((C+sig*U)'*W*Hsig(sig))))\...
    (Hsig(sig)'*W'*R(:,i) - Hsig(sig)'*(Wp1TC + sig*Wp1TU)*((CTU(sig))\((C+sig*U)'*R(:,i))));

X(:,i) = X(:,i) + W(:,1:m*s)*Y(:,i) + U*((CTU(sig)\((C+sig*U)'*(R(:,i) - W*Hsig(sig)*Y(:,i)))));

R(:,i) = R(:,i) - W*Hsig(sig)*Y(:,i) - (C+sig*U)*(((CTU(sig))\ (C+sig*U)'*(R(:,i) - W*Hsig(sig)*Y(:,i))));
end



iter = iter + 1;
resid(iter) = norm(R)

shift = shifts(shift_monitor);
[U] = block_harm_ritz(H,W,m,s,k,U,C,shift);
C = A*U;
[C,RR] = qr(C,0);
U = U/RR;

 if shift_recycle_method == 1
    shift_monitor = shift_monitor + 1;
    if shift_monitor == size(shifts,2) 
        shift_monitor = 1;
    end
 end


end