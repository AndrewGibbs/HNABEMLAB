function [X,Y] = rearrange( x,y )
%function designed to make tensor products of vectors x and y, and compute
%the corresponding weights W.
    Nx=length(x);   Ny=length(y);
    X=repmat(x,Ny,1);
    Y=reshape(repmat(y,1,Nx).',Ny*Nx,1);
end