function [z, Id0flag]=AllZeros(f,xmin,xmax,N)
% Inputs :
% f : function of one variable
% [xmin - xmax] : range where f is continuous containing zeros
% N : control of the minimum distance (xmax-xmin)/N between two zeros

Id0flag=false;

if (nargin<4)
    N=100;
end

thresh=1E-12;
XzeroTest=linspace(xmin,xmax,N);
if max(abs(f(XzeroTest)))<thresh
    Id0flag=true;
    z=[];
    return;
end

dx=(xmax-xmin)/N;
x2=xmin;
y2=f(x2);
z=[];
for i=1:N
    x1=x2;
    y1=y2;
    x2=xmin+i*dx;
    y2=f(x2);
    if (y1*y2<=0)                              % Rolle's theorem : one zeros (or more) present
        z=[z,fsolve(f,(x2*y1-x1*y2)/(y1-y2), optimoptions('fsolve','Display','off'))]; 
        % Linear approximation to guess the initial value in the [x1,x2] range.
    end
end