function x = quad_gauss_asymptotic(n)
%function x = quad_gauss_asymptotic(n)
%
%

x = zeros(floor(n/2),1);

for k=0:floor(n/2)-1
    theta_k = (4*(floor(n/2)-k)-1)/(4*n+2)*pi;
    x(k+1) = (1 - (n-1)/(8*n^3)-1/(384*n^4)*(39-28/sin(theta_k)^2)) * cos(theta_k);
end
x = sort(x);
if mod(n,2)==0
    x = [-flipud(x); x];
else
    x = [-flipud(x); 0; x];
end
