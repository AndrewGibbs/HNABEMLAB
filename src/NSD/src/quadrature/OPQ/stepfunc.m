% STEPFUNC Step function.
%
%    Given N and the Nx2 array ta of the knots and coefficients of the step function
%    s_N(x), this routine evaluates the step function at x, which may be vector-valued.
%
function y=stepfunc(N,x,ta)
if size(ta,1)~=N, error('array ta incorrectly sized'), end
y=zeros(size(x));
for n=1:N
  y=y+ta(n,2)*heaviside(ta(n,1)-x);
end
