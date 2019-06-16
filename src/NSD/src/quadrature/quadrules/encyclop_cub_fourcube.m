function [x,w] = encyclop_cub_fourcube(d, N, v)
%function [x,w] = encyclop_cub_fourcube(d, N, v)
%
%   Return a cubature formula of degree d and with N points on the four
%   cube. If there are more possibilities, v=2,3,... will select them.
%   Source: Encyclopedia of Cubature Formulae (R. Cools)

if nargin == 2
    v = 1;
end

switch d
    case 1
        switch N
            case 1
                x = gen_midpoint();
                w = 16;
        end
    case 3
        switch N
            case 8
                switch v
                    case 1
                        x = gen_centrallysymmetric_a(sqrt(3)/3);
                        w = ones(8,1) * 2;
                    case 2
                        x = gen_fullysymmetric_a000(1.15470053837925152901829756100391);
                        w = ones(8,1) * 2;
                end
            case 9
                x1 = gen_fullysymmetric_a000(1);
                w1 = ones(8,1) * 2.66666666666666666666666666666666;
                x2 = gen_midpoint();
                w2 = -5.33333333333333333333333333333333;
                x = [x1; x2];
                w = [w1; w2];
            case 16
                % This is actually a tensor-product rule...
                x = gen_fullysymmetric_a(0.577350269189625764509148780501957);
                w = ones(16,1);
        end
    case 5
        switch N
            case 24
                 x1 = gen_fullysymmetric_a000(0.894427190999915878563669467492510);
                 w1 = ones(8,1) * 1.11111111111111111111111111111111;
                 x2 = gen_fullysymmetric_a(0.707106781186547524400844362104849);
                 w2 = ones(16,1) * 0.444444444444444444444444444444444;
                 x = [x1; x2];
                 w = [w1; w2];
        end
    case 7
        switch N
            case 57
                x1 = [0 0 0 0];
                w1 = 0.613130128211648798876393271153710;
                x2 = gen_fullysymmetric_a000(0.986614654720553053315916432646997);
                w2 = ones(8,1) * 0.403855017620820896152581351450957;
                x3 = gen_fullysymmetric_aaa0(0.849521332456696146887976941995510);
                w3 = ones(32,1) * 0.157655993126334106135659249365851;
                x4 = gen_fullysymmetric_a(0.505408105853053508532167527628024);
                w4 = ones(16,1) * 0.444439871923693289722616246095711;
                x = [x1; x2; x3; x4];
                w = [w1; w2; w3; w4];
        end

end

end


function z = gen_midpoint()
    z = [0 0 0 0];
end

function z = gen_centrallysymmetric_a(a)
    z = [ a  a  a  a; -a -a -a -a;
          a  a -a -a; -a -a  a  a;
         -a  a -a  a;  a -a  a -a;
          a -a -a  a; -a  a  a -a];
end


function z = gen_fullysymmetric_a(a)
    z = [a a a a; a a a -a; a a -a a; a a -a -a; a -a a a; a -a a -a; a -a -a a; a -a -a -a; 
        -a a a a; -a a a -a; -a a -a a; -a a -a -a; -a -a a a; -a -a a -a; -a -a -a a; -a -a -a -a];
end

function z = gen_fullysymmetric_aaa0(a)
    z = [ a  a  a 0;  a  a  0  a;  a  0  a  a;  0  a  a  a;
          a  a -a 0;  a  a  0 -a;  a  0  a -a;  0  a  a -a;
          a -a  a 0;  a -a  0  a;  a  0 -a  a;  0  a -a  a;
          a -a -a 0;  a -a  0 -a;  a  0 -a -a;  0  a -a -a;
         -a  a  a 0; -a  a  0  a; -a  0  a  a;  0 -a  a  a;
         -a  a -a 0; -a  a  0 -a; -a  0  a -a;  0 -a  a -a;
         -a -a  a 0; -a -a  0  a; -a  0 -a  a;  0 -a -a  a;
         -a -a -a 0; -a -a  0 -a; -a  0 -a -a;  0 -a -a -a];
end

function z = gen_fullysymmetric_a000(a)
    z = [a 0 0 0; 0  a 0 0; 0 0  a 0; 0 0 0  a; ...
        -a 0 0 0; 0 -a 0 0; 0 0 -a 0; 0 0 0 -a];
end




