function z = apply_cub(f, x, w)
%function z = apply_cub(f, x, w)
%
%   Apply the cubature rule to f.

N = size(x,1);
z = 0;
for i = 1:N
    z = z + w(i) * f(x(i,1),x(i,2));
end
