function [x, y, r, w, bmx, ymb] = Graded_edges( a, b_m, b, b_p, c, angle, Sp, delta, N, x_width, y_width)
%Quadrature for functions supported near the same corner on separate sides.
%x(s),s\in[a,b_m]=S,    y(t),t\in[b_p,c]=T
%Sp=par.Sp; delta=par.Grad_delta; N=par.qppw;

if nargin==9
   x_width=b_m-a; y_width=c-b_p; 
end

    lamda_ab=(b-a)/(x_width); lambda_bc=(c-b)/(y_width);
    lamda=max(lamda_ab,lambda_bc);
    if abs(lamda-1)>1e-16 %account for rounding errors which can balls everything up and set p=-inf
        p_max= ceil(log_base( (lamda-1)/lamda, delta ));
        if ~isinteger(p_max)
            p_max=Sp;
        end
    else
        p_max=Sp;
    end
    
    %just define on [0,1]x[0,1] for now, can adjust after r has been
    %calculated.
    
    X=[];   W=[];
    [Z,w]=gauleg(N);
    Z=(Z+1)/2; w=w/2;
    X=delta^p_max*Z; W=delta^p_max*w;
    for p=fliplr(0:p_max-1)
        X=[X; delta^(p).*(delta+(1-delta)*Z)];
        W=[W; w*(1-delta)*delta^p];
    end
           
    [X,Y] = rearrange( X,X ); 
    [w1,w2] =rearrange(W,W);W=w1.*w2;

    %need to strectch, transpose, etc.
    bmx=(b-b_m)+X*(x_width);       ymb=(b_p-b)+Y*(y_width);
    x = b_m-X*(x_width);           y = b_p+Y*(y_width);
    w = W*(x_width)*(y_width);
    
    if angle~=pi
        r = sqrt(((b-b_m)+(x_width)*X).^2 + ((b_p-b)+(y_width)*Y).^2 ...
            -2*cos(angle)*((b-b_m)+(x_width)*X).*((b_p-b)+(y_width)*Y));
    else
        r=(b_p-b_m)+(x_width)*X + (y_width)*Y;
    end
end

%This function is one of four brothers. Corrections to one might be
%applicable to all.
    % Boxes_near
    % Boxes_near_EndArc
    % Boxes_touch
    % Boxes_touch_EndArc
    