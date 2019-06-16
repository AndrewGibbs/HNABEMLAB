% FEX31  A function used in Example 3.31.
%
function y=fex31(x,om)
n=size(x,1);
y=ones(n,1);
for k=1:n
  if x(k)~=0, y(k)=(pi*x(k)/om)/sin(pi*x(k)/om); end
end
