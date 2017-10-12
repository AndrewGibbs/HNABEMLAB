function [X,Y] = kron2( x,y )
%constructs repeated vectors such that each element of x and y match at a
%certain entry
    X=repmat(x(:),length(y),1);
    Y=repmat(y(:),length(x),1);
end

