function [resid] = unproj_rsbFOM(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method)


Y = zeros(m*s,s);
E = zeros(m*s,s);
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

[W,H,S] = block_Arnoldi(A,R,m,s,n);
E(1:s,1:s) = S;
WTU = W(:,1:m*s)'*U;
WTC = W(:,1:m*s)'*C;
UTU = U'*U;
UTC = U'*C;
UTWp1 = U'*W;


Hsig = @(zx) (zx*eyebar + H);

for i=1:s
sig = shifts(i);

Y(:,i) = ((sig*speye(m*s) + H(1:m*s,1:m*s)) - (WTC+sig*WTU)*(((UTC+sig*UTU)\(UTWp1*Hsig(sig)))))\...
    (E(:,i) - (WTC + sig*WTU)*(((UTC + sig*UTU)\(U'*R(:,i)))));

X(:,i) = X(:,i) + W(:,1:m*s)*Y(:,i) + U*((UTC+sig*UTU)\(U'*R(:,i) - UTWp1*Hsig(sig)*Y(:,i)));

R(:,i) = R(:,i) - W*Hsig(sig)*Y(:,i) - (C+sig*U)*( ((UTC + sig*UTU)\(U'*R(:,i) - UTWp1*Hsig(sig)*Y(:,i))));
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


end