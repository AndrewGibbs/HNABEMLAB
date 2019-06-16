function [c1,c2] = hoi_expf_F_sp1(s, a, fa, ga)
%function [c1,c2] = hoi_expf_F_sp1(s, a, fa, ga)
%
%   Return the full asymptotic expansion of F at a regular point a.
%   The expansion of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g.

% $Id$


c1 = zeros(s,1);
d1 = zeros(s,1);
c2 = zeros(s,1);
d2 = zeros(s,1);

I = 1i;

if s >= 1
    d1(1) = I*sqrt(2)/sqrt(I*ga(3));
    d2(1) = -I*sqrt(2)/sqrt(I*ga(3));
    c1(1) = -1/2*I*fa(1)*sqrt(2)*sqrt(pi)/sqrt(I*ga(3));
    c2(1) = -1/2*I*fa(1)*sqrt(2)*sqrt(pi)/sqrt(I*ga(3));
end

if s >= 2
    d(2) = -1/6*ga(4)*d(1)^2/ga(3);
    c(2) = 1/6*fa(1)*ga(4)*d(1)^2/ga(3)-I*fa(2)/ga(3);
end

% if s >= 3
%     d(3) = 
%     c(3) = 
% end
% 
% if s >= 4
%     d(4) = 
%     c(4) = 
% end
% 
% if s >= 5
%     d(5) = 
%     c(5) = 
% end
% 
% if s >= 6
%     d(6) = 
%     c(6) = 
% end
% 
% if s >= 7
%     d(7) = 
%     c(7) = 
% end
% 
% if s >= 8
%     d(8) = 
%     c(8) = 
% end
% 
% if s >= 9
%     d(9) = 
%     c(9) = 
% end
% 
% if s >= 10
%     d(10) = 
%     c(10) = 
% end
% 
% 
