% QUADGP General-purpose quadrature routine.
%
%    This is the in-house quadrature routine for use in mcdis and 
%    mccheb.
%
function xw=quadgp(N,i)
global mc uv AB
if(i>1&i<mc), xw=tr(uv,i); return, end
if mc==1
  if(AB(i,1)~=-Inf)&(AB(i,2)~=Inf)
    xw=tr(uv,i); return
  end
  if AB(i,1)~=-Inf, xw=rtr(uv,i); return, end
  if AB(i,2)~=Inf, xw=ltr(uv,i); return, end
  xw=str(uv,i); return
else
  if((i==1&AB(i,1)~=-Inf)|(i==mc&AB(i,2)~=Inf))
    xw=tr(uv,i); return
  end
  if i==1, xw=ltr(uv,i); return, end
end
xw=rtr(uv,mc);

function s=tr(t,i)
global AB 
s(:,1)=((AB(i,2)-AB(i,1))*t(:,1)+AB(i,2)+AB(i,1))/2;
s(:,2)=(AB(i,2)-AB(i,1)).*t(:,2).*wf(s(:,1),i)/2;

function s=str(t,i)
s(:,1)=t(:,1)./(1-t(:,1).^2);
s(:,2)=t(:,2).*wf(s(:,1),i).*(1+t(:,1).^2)./...
  ((1-t(:,1).^2).^2);

function s=rtr(t,i)
global AB 
s(:,1)=AB(i,1)+(1+t(:,1))./(1-t(:,1));
s(:,2)=2*t(:,2).*wf(s(:,1),i)./((1-t(:,1)).^2);

function s=ltr(t,i)
global AB 
s(:,1)=AB(i,2)-(1-t(:,1))./(1+t(:,1));
s(:,2)=2*t(:,2).*wf(s(:,1),i)./((1+t(:,1)).^2);

function y=wf(x,i)
y=exp(-x.*x);
