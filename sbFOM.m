function [resid] = sbFOM(A,B,X,shifts,m,s,n,tol)


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



while(norm(R) > tol)

[W,H,S] = block_Arnoldi(A,R,m,s,n);

E(1:s,1:s) = S;

Hsig = @(zx) H + zx*eyebar;
for i=1:s
sig = shifts(i);

Y(:,i) = (H(1:m*s,1:m*s) + sig*speye(m*s))\E(:,i);

X(:,i) = X(:,i) + W(:,1:m*s)*Y(:,i);

R(:,i) = R(:,i) - W*Hsig(sig)*Y(:,i);
end


iter = iter + 1;
resid(iter) = norm(R);

end


end