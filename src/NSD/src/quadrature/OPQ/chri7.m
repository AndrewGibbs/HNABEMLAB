% CHRI7  Modification by a special quadratic factor.
%
%    Given a weight function w(t) through the first n+1 recurrence
%    coefficients ab0 of its orthogonal polynomials, ab=CHRI7(n,
%    ab0,x) generates the first n (>1) recurrence coefficients of
%    the orthogonal polynomials relative to the modified weight 
%    function (t-x)^2 w(t). The alpha- and beta-coefficients of
%    the original weight function are to be provided in the first 
%    and second column of the (n+1)x2 array ab0; those of the 
%    modified weight function are returned in the first and second
%    column of the nx2 array ab.
%
function ab=chri7(N,ab0,x)
if N<2, error('N out of range'), end
N0=size(ab0,1); if N0<N+1, error('input array ab0 too short'), end
u=0; c=1; c0=0; s=0;
for k=1:N
  g=ab0(k,1)-x-u; cm1=c0; c0=c;
  if abs(c0)>5*eps
    p2=g^2/c0;
  else
    p2=cm1*ab0(k,2);
  end
  if k>1, ab(k,2)=s*(p2+ab0(k+1,2)); end
  s=ab0(k+1,2)/(p2+ab0(k+1,2));
  c=p2/(p2+ab0(k+1,2)); u=s*(g+ab0(k+1,1)-x);
  ab(k,1)=g+u+x;
end
ab(1,2)=ab0(1,2)*(ab0(2,2)+(x-ab0(1,1))^2);
