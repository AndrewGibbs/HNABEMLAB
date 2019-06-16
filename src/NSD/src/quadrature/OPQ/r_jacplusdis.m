% R_JACPLUSDIS Recurrence coefficients for polynomials orthogonal with
% respect to a normalized Jacobi weight function plus a discrete
% measure.
%
function ab=r_jacplus(N,alpha,beta,ty)
global mc mp iq idelta irout DM AB
global a b
a=alpha; b=beta;
mc=1; mp=size(ty,1); iq=1; idelta=2; Mmax=N+1;
DM=ty; AB=[-1 1]; eps0=1e3*eps;
[ab,Mcap,kount]=mcdis(N,eps0,@quadjp,Mmax);
