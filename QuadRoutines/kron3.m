function [X,Y,Z] = kron3( x,y,z )
%constructs repeated vectors such that each element of x and y match at a
%certain entry
    X=repmat(x(:),length(y)*length(z),1);
    Y=repmat(y(:),length(x)*length(z),1);
    Z=repmat(z(:),length(y)*length(x),1);
end