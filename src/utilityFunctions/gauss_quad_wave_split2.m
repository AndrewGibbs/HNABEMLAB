function [ x, w ] = gauss_quad_wave_split2(a, b, qppw, k,  ab_width )
    if nargin==4
        ab_width=abs(b-a);
    end
    
    wavelength = 2*pi/k;
    
    wavelengthsInWidth = ceil(ab_width/wavelength);
    
    sub_width = ab_width/wavelengthsInWidth;
    
    x=[]; w=[];
    
    for q=1:wavelengthsInWidth
        [x_,w_]=gauss_quad((a+(q-1)*sub_width),(b-(wavelengthsInWidth-q)*sub_width),qppw);
        x=[x; flipud(x_)]; w=[w; flipud(w_)];
    end
end