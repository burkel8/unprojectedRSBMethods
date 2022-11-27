%Function to compute approximation to f(A)b when f is the inverse square root function
% using implementation 2 of r(FOM)^2. The quadrature rule used is outlined
% in section 7 of the preprint.

%Inputs: b - vector b for which we want to approximate f(A)b
%        V,H - basis of Krylov subspace and hessenberg matrix H constructed
%        from Arnoldi process to build basis for Km(A,b)
%        m - dimension of Krylov subspace
%        k - dimension of recycling subspace
%        U - recycling subspace
%        C - matrix C = A*U.
%        num_quad - number of quadrature points to be used

%Output: deflated_approx an n x 1 vector storing the approximation to f(A)b
function [deflated_approx] = rFOM2_v2_invSqrt(b,V,H,m,k,U,C,num_quad)

term1 = zeros(m+k,1);

 % Define constant factors appearing in quadrature integration
 [U,D] = scale_cols_of_U(U,k);
 Vhat = [U V(:,1:m)];
 What = [C V(:,1:m)];
 G = zeros(m+k,m+k);
 G(1:k,1:k) = D;
 G(k+1:m+k,k+1:m+k) = H(1:m,1:m);
 UmC = U-C;
 e = zeros(m,1);
 e(m)=1;
 hterm = H(m+1,m)*V(:,m+1)*e';
 VTb = Vhat'*b;
 VTW = Vhat'*What;

 R = @(zx) [zx*UmC  -hterm];
yy = @(zx) (VTW*(zx*speye(m+k)-G) + Vhat'*R(zx))\VTb;

const = -2/pi;

%compute quadrature nodes and weights
weights = pi/num_quad*ones(1,num_quad);
t = zeros(1,num_quad);
for ii = 1:num_quad
t(ii) = cos((2*ii-1)/(2*num_quad) * pi);
end
tt = -1*(1-t)./(1+t);

%perform quadrature
for j = 1:num_quad
  yterm = yy(tt(j));
  term1 = term1 + weights(j)*(1/(1+t(j)))*yterm;
end

%Compute approximation
deflated_approx = Vhat*const*term1;
end