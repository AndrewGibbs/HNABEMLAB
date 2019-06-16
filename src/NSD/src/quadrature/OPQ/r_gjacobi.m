% R_GJACOBI Recurrence coefficients for orthogonal polynomials
% relative to a generalized Jacobi weight function.
%
%    The generalized Jacobi weight function is specified by mc+1
%    nodes -1, t_2,...,t_mc,1 and corresponding exponents
%    b,g_2,...,g_mc,a. These are assumed available in the (mc+1)x2
%    array tc, the nodes occupying the first column of tc and
%    the exponents the second column. It is assumed that the
%    calling program selects between the Stieltjes and Lanczos
%    procedure by setting the global variable irout accordingly
%    and also specifies the error tolerance through the global
%    variable eps0. The call ab=R_GJACOBI(n,a) is the same as
%    ab=R_GJACOBI(n,a,a) and ab=R_GJACOBI(n) the same as
%    ab=R_GJACOBI(n,0,0).
%
function ab=r_gjacobi(N,a,b)
global mc mp iq idelta irout eps0 AB
global tc Mcap kount
if nargin<2, a=0; end; if nargin<3, b=a; end
mc=size(tc,1)-1; mp=0; iq=1; idelta=2;
Mmax=N+200; 
if mc<2, error('mc must be at least equal to 2'), end
AB=[[tc(1:mc,1) tc(2:mc+1,1)]]; 
[ab,Mcap,kount]=mcdis(N,eps0,@quadgj,Mmax);
