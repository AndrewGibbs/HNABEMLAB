function [ x,y,w] = Gauss_2D_quad_widths( A,C, N, A_width, C_width)
%Using Steve's notes to create a general 2D quadrature

    if nargin==3
       A_width=A(2)-A(1);
       C_width=C(2)-C(1);
    end
    
    [x_,w_x]=gauss_quad(A(1),A(2),N);
    [y_,w_y]=gauss_quad(C(1),C(2),N);
        [x,y] = rearrange(x_, y_);
        [Wx,Wy] = rearrange( w_x, w_y );
        w=Wx.*Wy;
end