% R_KRONROD Jacobi-Kronrod matrix.
%
%    ab=R_KRONROD(n,ab0) produces the alpha- and beta-elements in 
%    the Jacobi-Kronrod matrix of order 2n+1 for the weight
%    function (or measure) w. The input data for the weight
%    function w are the recurrence coefficients of the associated
%    orthogonal polynomials, which are stored in the array ab0 of
%    dimension [ceil(3*n/2)+1]x2. The alpha-elements are stored in
%    the first column, the beta-elements in the second column, of 
%    the (2*n+1)x2 array ab.
%
%    Supplied by Dirk Laurie, 6-22.1998
%
function ab=r_kronrod(N,ab0)
if length(ab0)<ceil(3*N/2)+1, error('array ab0 too short'), end
a=zeros(2*N+1,1); b=a;
k=0:floor(3*N/2); a(k+1)=ab0(k+1,1);
k=0:ceil(3*N/2); b(k+1)=ab0(k+1,2);
s=zeros(floor(N/2)+2,1); t=s; t(2)=b(N+2);
for m=0:N-2,
  k=floor((m+1)/2):-1:0; l=m-k;
  s(k+2)=cumsum((a(k+N+2)-a(l+1)).*t(k+2)+b(k+N+2).*s(k+1)-b(l+1).*s(k+2));
  swap=s; s=t; t=swap;
end
j=floor(N/2):-1:0; s(j+2)=s(j+1);
for m=N-1:2*N-3,
  k=m+1-N:floor((m-1)/2); l=m-k; j=N-1-l;
  s(j+2)=cumsum(-(a(k+N+2)-a(l+1)).*t(j+2)-b(k+N+2).*s(j+2)+b(l+1).*s(j+3));
  j=j(length(j)); k=floor((m+1)/2);
  if rem(m,2)==0, a(k+N+2)=a(k+1)+(s(j+2)-b(k+N+2)*s(j+3))/t(j+3);
  else b(k+N+2)=s(j+2)/s(j+3);
  end
  swap=s; s=t; t=swap;
end
a(2*N+1)=a(N)-b(2*N+1)*s(2)/t(2);
ab=[a b];
