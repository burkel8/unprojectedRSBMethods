%Function to compute approximation to f(A)b when f is the inverse square root function
% using implementation 1 of r(FOM)^2. The quadrature rule used is outlined
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
function [deflated_approx] = rFOM2_v1_invSqrt(b,V,H,m,k,U,C,num_quad)

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
  yterm = y(tt(j));
  term1 = term1 + (1/(1+t(j)))*weights(1,j)*yterm;
  term2 = term2 + (1/(1+t(j)))*weights(1,j)*z1(tt(j),yterm);
end

%Compute approximation
deflated_approx = const*V(:,1:m)*term1 + const*U*term2;
end