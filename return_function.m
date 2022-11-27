%Function which returns the function to be used in experiments.
%User can add a function here if they wish.

function [f_scalar, f_matrix] = return_function(problem)

f1_scalar = @(zx) 1.0/zx;
f1_matrix = @(Ax,bx) Ax\bx;

f2_scalar = @(zx) 1./sqrt(zx);
f2_matrix = @(Ax,bx) sqrtm(full(Ax))\bx;

f3_scalar = @(zx) log(zx);
f3_matrix = @(Ax,bx) logm(full(Ax))*bx;

f4_scalar = @(zx) sqrt(zx);
f4_matrix = @(Ax,bx) sqrtm(full(Ax))*bx;

if problem =="inverse"
    f_scalar = @(zx) f1_scalar(zx);
    f_matrix = @(Ax,bx) f1_matrix(Ax,bx);
elseif problem == "invSqrt"
    f_scalar = @(zx) f2_scalar(zx);
    f_matrix = @(Ax,bx) f2_matrix(Ax,bx);
elseif problem == "log"
    f_scalar = @(zx) f3_scalar(zx);
    f_matrix = @(Ax,bx) f3_matrix(Ax,bx);
elseif problem == "sqrt"
    f_scalar = @(zx)  f4_scalar(zx);
    f_matrix = @(Ax,bx) f4_matrix(Ax,bx);
else
    error("ERROR : unknown function chosen!\n");
end


end