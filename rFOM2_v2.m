%Function to compute approximation to f(A)b using implementation 2 of
%r(FOM)^2. The quadrature rule used is a trapezoidal rule.

%Inputs: b - vector b for which we want to approximate f(A)b
%        V,H - basis of Krylov subspace and hessenberg matrix H constructed
%        from Arnoldi process to build basis for Km(A,b)
%        m - dimension of Krylov subspace
%        k - dimension of recycling subspace
%        U - recycling subspace
%        C - matrix C = A*U.
%        num_quad - number of quadrature points to be used
%        f_scalar - the scalar form of the matrix function f(z) , z scalar

%Output: deflated_approx an n x 1 vector storing the approximation to f(A)b
function [deflated_approx] = rFOM2_v2(b,V,H,m,k,U,C,num_quad, f_scalar)

%Construct circular contour with radius r and centre circle_centre such
%that the spectrum of A is contained within the contour.
s = eigs(H(1:m,1:m),1,'smallestreal');
l = eigs(H(1:m,1:m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

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

 R = @(zx) [zx*UmC  -hterm]; %matrix which depends on quadrature point

% function yy representing solution (6.4) on preprint.
yy = @(zx) (VTW*(zx*speye(m+k)-G) + Vhat'*R(zx))\VTb;

delta_theta = 2*pi / num_quad; %distance between quadrature points.
const = r/num_quad;

%perform quadrature
for j = 1:num_quad

  theta = (j-1)*delta_theta;

  %move to new quadrature point.
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = yy(z);
  term1 = term1 + common_factor*yterm;
end

%Compute approximation
deflated_approx = Vhat*const*term1;
end