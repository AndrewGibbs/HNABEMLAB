% EPSALG Epsilon algorithm.
%
%    Given a (column) vector s of dimension n, this routine applies
%    to it the epsilon algorithm, returning the results in the
%    nx(n+1) array E. The values of interest are located in the
%    even-numbered columns of E. Preferably, n is to be chosen
%    an odd integer.
%
function E=epsalg(n,s)
if n<2, error('n too small'), end
if size(s,1)~=n, error('s and n incompatible'), end
E=zeros(n,n+1);
E(:,2)=s;
for k=3:n+1
  for m=1:n+2-k
    if E(m+1,k-1)==E(m,k-1), break, end
    E(m,k)=E(m+1,k-2)+1/(E(m+1,k-1)-E(m,k-1));
  end
end
