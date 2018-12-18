function [ x, w ] = gauss_quad_wave_split(a, b, qppw, k,  ab_width )
    if nargin==4
        ab_width=abs(b-a);
    end
    
    wavelength = 2*pi/k;
    
    wavelengthsInWidth = ceil(width/wavelenth);
    
    sub_width = ab_width/wavelengthsInWidth;
    
    x=[]; w=[];
    
    for q=1:subs
        [x_,w_]=gauss_quad((a+(q-1)*sub_width),(b-(subs-q)*sub_width),qppw);
        x=[x.' fliplr(x_.')].'; w=[w.' fliplr(w_.')].';
    end
end