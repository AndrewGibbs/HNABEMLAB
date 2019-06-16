function [x,w] = encyclop_cub_triangle(d, N, v)
%function [x,w] = encyclop_cub_triangle(d, N, v)
%
%   Return a cubature formula of degree d and with N points on the unit
%   simplex in 2d. If there are more possibilities, v=2,3,... will select
%   them.
%   Source: Encyclopedia of Cubature Formulae (R. Cools)

if nargin == 2
    v = 1;
end

switch d
    case 1
        switch N
            case 1
                x = gen_midpoint();
                w = 1/2;
        end
    case 2
        switch N
            case 3
                x = gen_fullysymmetric_aab(1/6, 0);
                w = ones(3,1) * 1/6;
        end
    case 3
        switch N
            case 4
                x1 = gen_midpoint();
                w1 = -0.28125;
                x2 = gen_fullysymmetric_aab(0.2, 0);
                w2 = ones(3,1) * 0.260416666666666666666666666666666;
                x = [x1; x2];
                w = [w1; w2];
        end
    case 5
        switch N
            case 7
                x1 = gen_midpoint();
                w1 = 0.1125;
                x2 = gen_fullysymmetric_aab(0.101286507323456338800987361915123, 0);
                w2 = ones(3,1) * 0.0629695902724135762978419727500906;
                x3 = gen_fullysymmetric_aab(0.470142064105115089770441209513447, 0);
                w3 = ones(3,1) * 0.0661970763942530903688246939165759;
                x = [x1; x2; x3];
                w = [w1; w2; w3];
        end
end

end

function z = gen_midpoint()
    v1 = [0 0];
    v2 = [1 0];
    v3 = [0 1];
    z = 1/3*(v1+v2+v3);
end


function z = gen_fullysymmetric_aab(a, b)
    v1 = [0 0];
    v2 = [1 0];
    v3 = [0 1];
    b = 1-2*a;
%     a = a/c; b = b/c;
    z = [a*v1+a*v2+b*v3; a*v1+b*v2+a*v3; b*v1+a*v2+a*v3];
end

