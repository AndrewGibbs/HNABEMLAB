function [x,w] = LinearNSD(Twidth, kwave, alphaPhase, xiSing, Npts)
    %parameters which determine cases, determined experimentally:
    c1 = 2*pi;
    c2 = .2;
    
    c3 = c1*c2; % local region (scaled by frequency) which contributes majority of oscillatory integral
    
    %now categorise integral:
    if alphaPhase*kwave*Twidth<= c1 %(non-oscillatory)
        oscillatory = false;
        if abs(xiSing)/Twidth >= c2 %(non-singular)
            singular = false;
        else
            singular = true; %(singular)
        end
    else %(oscillatory)
        oscillatory = true;
       if abs(xiSing)*(kwave*alphaPhase) >= c3 %(non-singular)
           singular = false;
       else
           singular = true; %(singular)
       end
    end
    
    %now split into four cases:
    if oscillatory %oscillatory
        if singular %singular
            [x0,w0] = LaplaceSD(xiSing, kwave*alphaPhase, 2*Npts, true);%GGLa
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
