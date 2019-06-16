% MOLL Mollifier function.
%
function m=moll(t)
m=abs(t);
m(find(m<1))=1;
