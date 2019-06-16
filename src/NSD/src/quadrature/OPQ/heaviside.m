% HEAVISIDE The Heaviside function.
%
function y=heaviside(x)
y=ones(size(x)); y(find(x<0))=0;
