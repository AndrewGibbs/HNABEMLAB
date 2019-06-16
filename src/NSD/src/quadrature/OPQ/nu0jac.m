% NU0JAC Initial backward recurrence index for Jacobi weight.
%
%    The call nu0=NU0JAC(n,z,eps0) computes the initial index nu0
%    in the continued fraction algorithm for generating to a
%    relative accuracy eps0 the Cauchy integral int_{-1}^1
%    pi_n(t)w(t)dt/(z-t) for large n, and complex z outside of
%    [-1,1]. Here, pi_n is the Jacobi polynomial of degree n, and
%    w the corresponding weight function on [-1,1].
%
%    See also CAUCHY and KERNEL.
%
function nu0=nu0jac(N,z,eps0)
num=log(1/eps0);
den=2*(log(abs(z+(z-1).^(1/2).*(z+1).^(1/2))));
nu0=fix(N+1+num./den);
