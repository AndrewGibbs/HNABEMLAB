function [ x,w ] = pathQuadFinite( b,S, m, freq,N )
%special quadrature rule for finite paths, basically just a Gaussian
%quadrature rule, but this might change

if ~isempty(S)
    error('Cant handle singularities on finite paths yet');
end
    
        [ x, w_preExp ] =  quad_gauss(N, 0, b);
        %no rescaling of this integral:
        w=w_preExp.*exp(-freq*x.^m);

end

