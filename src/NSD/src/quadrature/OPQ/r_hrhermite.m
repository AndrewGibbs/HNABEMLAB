% R_HRHERMITE  Half-range Hermite measure.
%
%    The call ab=R_HRHERMITE(N,Mmax,eps0) generates the first
%    N recurrence coefficients for the half-range Hermite measure,
%    the alpha- and beta-coefficients being stored in the first
%    resp. second column of the Nx2 array ab. The routine uses
%    the multi-component discretization routine mcdis.m, which
%    requires the parameter Mmax to control the discretization
%    parameter M, and the parameter eps0 to control the accuracy.
%
function ab=r_hrhermite(N,Mmax,eps0)
global mc mp iq idelta irout AB
mc=4; mp=0; iq=0; idelta=1; irout=1;
AB=[[0 3];[3 6];[6 9];[9 Inf]];
[ab,Mcap,kount]=mcdis(N,eps0,@quadgp,Mmax);
