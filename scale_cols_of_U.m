%A function which takes in a U and returns Unew and D such Unew = UD has
%unit columns.

%Input: U - Tall Skinny matrix we wish to scale.
%       k - Number of columns of U

%Output U such that U = UD;
%       D - k x k scaling matrix 

function [U,D] = scale_cols_of_U(U,k)
D = zeros(k,k);
for i = 1:k
   D(i,i) = 1.0/norm(U(:,i));
   U(:,i) = U(:,i)*D(i,i);
end
end