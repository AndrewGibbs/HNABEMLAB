% NU0HER Initial backward recurrence index for Hermite weight.
%
%    The call nu0=NU0HER(n,z,eps0) computes the initial index nu0
%    in the continued fraction algorithm for generating to a
%    relative accuracy eps0 the Cauchy integral int_{-Inf}^Inf
%    pi_n(t)w(t)dt/(z-t) for large n, and complex z outside of
%    [-Inf,Inf]. Here, pi_n is the Hermite polynomial of degree n,
%    and w the correponding weight function.
%
%    See also CAUCHY and KERNEL.
%
function nu0=nu0her(N,z,eps0)
nu0=fix(2*(sqrt((N+1)/2)+.25*log(1/eps0)./abs(imag(z))).^2);
