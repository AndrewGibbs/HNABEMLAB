function [x,w] = LinearNSD(Twidth, kwave, alphaPhase, xiSing, Npts)
    %parameters which determine cases, determined experimentally:
    c_osc = 2;
    c_sing = .5;

    if kwave*alphaPhase*Twidth/(2*pi)<=c_osc
        oscillatory = false;
        if abs(xiSing)/Twidth >= c_sing
            singular = false;
        else
            singular = true;
        end
    else
        oscillatory = true;
        if abs(xiSing)*kwave*alphaPhase/(2*pi*c_osc) >= c_sing
            singular = false;
        else
            singular = true;
        end
    end
    
    %now split into four cases:
    if oscillatory %oscillatory
        if singular %singular
            [x0,w0] = LaplaceSD(xiSing, kwave*alphaPhase, Npts, true);%GGLa
            [x2,w2] = LaplaceSD(Twidth, kwave*alphaPhase, Npts, false); %GLa
            if xiSing == 0
                x = [x0; x2]; 
                w = [w0; -w2];
            else
                [x1, w1] = genGaussLog( 2*Npts, xiSing, 0); %-GGLeg
                x = [x0; x1; x2];
                w = [w0; -w1.*exp(1i*kwave*alphaPhase*x1); -w2];
            end
        else %not singular
            [x0,w0] = LaplaceSD(0, kwave*alphaPhase, Npts, false);%GLa
            [x1,w1] = LaplaceSD(Twidth, kwave*alphaPhase, Npts, false);%GLa
            x = [x0; x1];
            w = [w0; -w1];
        end
    else %not oscillatory
        if singular
            [x1,w1] = genGaussLog( 2*Npts, xiSing, Twidth);%GGLeg
            if xiSing == 0
                x = x1;
                w = w1.*exp(1i*kwave*alphaPhase*x);
            else
                [x0,w0] = genGaussLog( 2*Npts, xiSing, 0);%-GGLeg
                x = [x0; x1];
                w = [-w0.*exp(1i*kwave*alphaPhase*x0); w1.*exp(1i*kwave*alphaPhase*x1)];
            end
        else %not singular
            [x_, w_] = quad_gauss(Npts, 0, Twidth); %GLeg
            x = x_;
            w = w_.*exp(1i*kwave*alphaPhase*x);
        end
    end
end
