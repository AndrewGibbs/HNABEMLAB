function [x,w] = quad_hp(a, b, n, sigma, mu, minsize)
%function [x,w] = quad_hp(a, b, n, sigma, mu, minsize)
%
%   Return the points and weights of a composite quadrature rule that has a
%   hp-type grading towards the point a.

% $Id$

if nargin == 5
    minsize = 0;
end

x = [];
w = [];
t2 = 1;
t1 = 1;
for i = n:-1:1
    t1 = t2 * sigma;
    p = max(2, floor(mu*n) + 1);
    u = t1*(b-a)+a;
    v = t2*(b-a)+a;
    if abs(u-v) < minsize
        % undo the last change and exit the loop
        t1 = t1/sigma;
        break
    end
    [x_new,w_new] = quad_gauss(p, u, v);
    x = [x; x_new];
    w = [w; w_new];
    t2 = t1;
end

[x_new,w_new] = quad_gauss(2, a, t1*(b-a)+a);
x = [x; x_new];
w = [w; w_new];
