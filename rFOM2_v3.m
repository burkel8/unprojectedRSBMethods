%Function to compute approximation to f(A)b using implementation 3 of
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
%        f_matrix - the matrix function f(A), for n x n matrix A.
%Output: deflated_approx an n x 1 vector storing the approximation to f(A)b
function [deflated_approx] = rFOM2_v3(b,V,H,m,k,U,C,num_quad, f_scalar, f_matrix)

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

em = zeros(m,1); em(m) = 1.0;
VTb = Vhat'*b;
VTW = Vhat'*What;
VTWinvVTb = (VTW)\VTb;

hterm = -H(m+1,m)*V(:,m+1)*em';
I = speye(m+k);

R = @(zx) [zx*UmC,hterm];
Gz = @(zx) VTW*(zx*speye(k+m)-G) ;
B = @(zx) Vhat'*R(zx);

%Function yy representing the quadrature integral I in eq (6.8) of preprint.
yy = @(zx) Gz(zx)\((I + B(zx)*(Gz(zx)\I))\(B(zx)*(Gz(zx)\VTb)));

delta_theta = 2*pi / num_quad;
const = r/num_quad;

%perform quadrature
for j = 1:num_quad

  theta = (j-1)*delta_theta;

  %move to new point in circle.
  z = r*exp(1i*theta) + circle_centre;
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = yy(z);
  term1 = term1 + common_factor*yterm;
end

%Compute Approximation. First term is closed form term in eq (6.8) in
%preprint. Second term comes from the quadrature.
deflated_approx = Vhat*f_matrix(G,VTWinvVTb) - const*Vhat*term1;
end