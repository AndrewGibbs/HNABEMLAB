% NU0LAG Initial backward recurrence index for Laguerre weight.
%
%    The call nu0=NU0LAG(n,z,a,eps0) computes the initial index 
%    nu0 in the continued fraction algorithm for generating to a
%    relative accuracy eps0 the Cauchy integral int_0^Inf
%    pi_n(t)w(t)dt/(z-t) for large n, and complex z outside of
%    [0,Inf]. Here, pi_n is the generalized Laguerre  polynomial
%    of degree n with parameter a, and w the corresponding 
%    Laguerre weight function.
%
%    See also CAUCHY and KERNEL.
%
function nu0=nu0lag(N,z,a,eps0)
num=log(1/eps0);
den=4*imag(sqrt(z));
nu0=fix((sqrt(N+(a+3)/2)+num./den).^2-(a+1)/2);
