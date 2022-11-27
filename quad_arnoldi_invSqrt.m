%This function computes the standard Arnoldi approximation using for f(A)b
%where f is the inverse square root function. A brief outline of the
%quadrature rule is given in section 7 of the preprint and references
%therin.

%Inputs: 
%        V,H - basis of Krylov subspace and hessenberg matrix H constructed
%        from Arnoldi process to build basis for Km(A,b)
%        m - dimension of Krylov subspace
%        num_quad - number of quadrature points to be used

%Output: standard_approx - an n x 1 vector storing the aproximation to f(A)b

function [standard_approx] = quad_arnoldi_invSqrt(V,H,m,num_quad)

e1 = zeros(m,1);
e1(1) = 1;
term1 = zeros(m,1);

%Function y representing the FOM approximation to the shifted linear system
%(sig*I - A)x(sig) = b
y = @(zx) (zx*speye(m)-H(1:m,1:m))\e1;

const = -2/pi;

%Compute quadrature nodes and weights
weights = pi/num_quad*ones(1,num_quad);
t = zeros(1,num_quad);
for ii = 1:num_quad
t(ii) = cos((2*ii-1)/(2*num_quad) * pi);
end
tt = -1*(1-t)./(1+t);

%perform quadrature
for j = 1:num_quad
   term1 = term1 + weights(j)*(1/(1+t(j)))*y(tt(j));
end

%Compute approximation
standard_approx = const*V(:,1:m)*term1;
end