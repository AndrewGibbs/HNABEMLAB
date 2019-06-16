function [ x, w ] = gauss_quad_wave_split(a, b, subs, N, ab_width )
    if nargin==4
        ab_width=abs(b-a);
    end

    %subs=ceil(2*ab_width/wavelength);
    sub_width=ab_width/subs;
    
    x=[]; w=[];
    
    for q=1:subs
        [x_,w_]=gauss_quad((a+(q-1)*sub_width),(b-(subs-q)*sub_width),qppw);
        x=[x.' fliplr(x_.')].'; w=[w.' fliplr(w_.')].';
    end
end