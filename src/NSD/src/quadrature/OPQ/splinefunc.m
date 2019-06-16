% SPLINEFUNC Spline function.
%
%    Given N, M, and the Nx2 array ta of the knots and coefficients of the spline,
%    this routine evaluates the spline function s_{N,M}(x) at x, which may be
%    vector-valued.
%
function y=splinefunc(N,M,x,ta)
if M<1, error('M too small'), end
if size(ta,1)<N, error('array ta too short'), end
y=zeros(size(x));
for n=1:N
  y=y+ta(n,2)*((ta(n,1)-x).*heaviside(ta(n,1)-x)).^M;
end
