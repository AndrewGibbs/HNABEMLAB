clear all

syms x1 x2 yd1 yd2 P1 P2 d1 d2 pm cosTheta ty_dist sx_dist r L;

assume((yd1^4 + 2*yd1^2*yd2^2 - yd1^2 + yd2^4 - yd2^2) ==0);
assume(pm^2==1)
assume(L>0)
%assume(- pm^2*yd1^2 - pm^2*yd2^2 + yd1^4 + 2*yd1^2*yd2^2 + yd2^4 == 0);
%assume((2*(P1*yd1^3 + P2*yd2^3 - x1*yd1^3 - x2*yd2^3 + P1*yd1*yd2^2 + P2*yd1^2*yd2 + pm^2*x1*yd1 + pm^2*x2*yd2 - x1*yd1*yd2^2 - x2*yd1^2*yd2 - P1*pm^2*yd1 - P2*pm^2*yd2)) ==0 , - pm^2*yd1^2 - pm^2*yd2^2 + yd1^4 + 2*yd1^2*yd2^2 + yd2^4 == 0);

g1a = symfun(sqrt( (ty_dist + r)^2 + sx_dist^2  -2*(ty_dist + r)*sx_dist*cosTheta) + pm*(ty_dist + r),r);
g1b = symfun(sqrt( (ty_dist + r)^2 + sx_dist^2  -2*(ty_dist + r)*sx_dist*cosTheta) + pm*(L - (ty_dist + r)),r);
g2a = symfun(sqrt( (ty_dist + r)^2 + sx_dist^2  -2*(ty_dist + r)*sx_dist*cosTheta) + (P1*d1 + P2*d2 + (ty_dist + r)*(d1*yd1 + d2*yd2)),r);
g2b = symfun(sqrt( (ty_dist + r)^2 + sx_dist^2  -2*(ty_dist + r)*sx_dist*cosTheta) + (P1*d1 + P2*d2 + (L - (ty_dist + r))*(d1*yd1 + d2*yd2)),r);
dg1a = diff(g1a,r);
dg1b = diff(g1b,r);
dg2a = diff(g2a,r);
dg2b = diff(g2b,r);
xi1a = solve(dg1a == 0,r);
xi1b = solve(dg1b == 0,r);
xi2a = solve(dg2a == 0,r);
xi2b = solve(dg2b == 0,r);
assume(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2==1);
assume(d1*yd1 + d2*yd2==1)
dg2a = diff(g2a,r);
dg2b = diff(g2b,r);
xi2aGraze = solve(dg2a == 0,r);
xi2bGraze = solve(dg2b == 0,r);
% collect(dg3,r)
% 
% (2*(yd1^2 + yd2^2)*(yd1*(P1 - x1) + yd2*(P2 - x2)) - 2*yd2*(P2 - x2) - 2*yd1*(P1 - x1))*t + (yd1*(P1 - x1) + yd2*(P2 - x2))^2 - (P1 - x1)^2 - (P2 - x2)^2
% 
% r = symfun(sqrt((x1 - (P1 + t*yd1))^2+(x2 - (P2 + t*yd2))^2),t);
% solve(r==0,t)

%summary: somehow, the S_l\varphi_n(x_m) phase has no stationary points,
%but S_l\Psi_n(x_m) does (EXCEPT in the case of grazing). These are:

%for singularity close to a:
%   (d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*sx_dist + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + cosTheta*d1^2*sx_dist*yd1^2 + cosTheta*d2^2*sx_dist*yd2^2 + 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1) - ty_dist
%   - ty_dist - (cosTheta*sx_dist + d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*d1^2*sx_dist*yd1^2 - cosTheta*d2^2*sx_dist*yd2^2 - 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1)
%and if close to b:
%    (d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*sx_dist + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + cosTheta*d1^2*sx_dist*yd1^2 + cosTheta*d2^2*sx_dist*yd2^2 + 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1) - ty_dist
%  - ty_dist - (cosTheta*sx_dist + d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*d1^2*sx_dist*yd1^2 - cosTheta*d2^2*sx_dist*yd2^2 - 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1)
 