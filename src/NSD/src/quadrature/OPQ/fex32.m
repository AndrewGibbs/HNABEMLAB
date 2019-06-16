% FEX32  A function used in Example 3.32.
%
function y=fex32(x,om)
n=size(x,1);
y=zeros(n,1);
for k=1:n
  y(k)=gamma(1+x(k))/(x(k)+om);
end
