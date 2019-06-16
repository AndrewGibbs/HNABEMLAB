% R_MODBESS Recurrence coefficients for a weight function 
% involving a modified Bessel function.
%
%    This generates the first N recursion coefficients ab
%    for the weight function w(t)=t^a*K_0(t) on [0,Inf],
%    where K_0 is the modified Bessel function of order 0.
%    The input parameter Mmax is the maximum allowed M-value 
%    in the multiple-component discretization routine mcdis,
%    and eps0 an error tolerance also used by mcdis.
%
function ab=r_modbess(N,a,Mmax,eps0)
global mc mp iq idelta irout AB
global abjac abjaclog ablag
if (N<=0)|(a<=-1), error('parameter(s) out of range'), end
mc=3; mp=0; iq=1; idelta=2; irout=1;
AB=[[0 1];[0 1];[0 Inf]];
abjac=r_jacobi01(Mmax,0,a);
abjaclog=r_jaclog(Mmax,a);
ablag=r_laguerre(Mmax);
ab=mcdis(N,eps0,@quadbess,Mmax);
