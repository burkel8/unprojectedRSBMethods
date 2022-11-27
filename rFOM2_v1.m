%Function to compute approximation to f(A)b using implementation 1 of
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
function [deflated_approx] = rFOM2_v1(b,V,H,m,k,U,C,num_quad, f_scalar)

%Construct circular contour with radius r and centre circle_centre such
%that the spectrum of A is contained within the contour.
s = eigs(H(1:m,1:m),1,'smallestreal');
l = eigs(H(1:m,1:m),1,'largestreal');
sm_eig = s ;
lg_eig = l;
shift = 0;
circle_centre = sm_eig+(lg_eig-sm_eig)/2; % for f(z) = 1.0/z
r = (sm_eig/2)+(circle_centre-sm_eig) + shift;     

term1 = zeros(m,1);
term2 = zeros(k,1);

%Define constant factors appearing in quadrature integration
VTU = V(:,1:m)'*U;
VTC = V(:,1:m)'*C;
UTU = U'*U;
UTC = U'*C;
UTV = U'*V(:,1:m);
VTb = V(:,1:m)'*b;
UTb = U'*b;
UTVp1H = U'*V*H;

%function y representing solution to linear system (6.2) in preprint
y = @(zx) (zx*speye(m) - H(1:m,1:m) - (zx*VTU - VTC)*( (zx*UTU-UTC)\ ...
    (zx*UTV - UTVp1H) ))\(VTb - (zx*VTU - VTC)*( (zx*UTU - UTC)\UTb));

% function z1 representing the z correction in terms of y
z1 = @(zx,yx) (zx*UTU - UTC)\(UTb - (zx*UTV - UTVp1H)*yx) ;

delta_theta = 2*pi / num_quad;
const = r/num_quad;

%Perform quadrature
for j = 1:num_quad

  theta = (j-1)*delta_theta;

  %move to new quadrature point
  z = r*exp(1i*theta) + circle_centre;
  
  common_factor = f_scalar(z)*exp(1i*theta);
  yterm = y(z);
  term1 = term1 + common_factor*yterm;
  term2 = term2 + common_factor*z1(z,yterm);
end

%Compute approximation
deflated_approx = const*V(:,1:m)*term1 + const*U*term2;
end