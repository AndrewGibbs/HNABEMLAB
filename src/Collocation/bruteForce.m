function I = bruteForce(Op, fun, funSide, colPt)
% use matlab's standard quadrature routine for integrals we haven't tackled
% appropriately yet:
    qppw = 15;
    supp = fun.getSupp(funSide);
    a = supp(1);
    b = supp(2);
    F = @(y) Op.kernel(colPt.x, y, colPt.side, funSide) .* fun.eval(y);
    if colPt.side == funSide
        c = colPt.x;
        I_a = integral(F, a, b, 'arrayValued', true,'RelTol',1e-8);
        I_c = integral(F, b, c, 'arrayValued', true,'RelTol',1e-8);
        I = I_a + I_c;
    else
        colPt_R2 = Op.domain.component(colPt.side).trace(colPt.x);
        a_R2 = Op.domain.component(funSide).trace(a);
        %b_R2 = Op.domain.component(funSide).trace(b);
        %get parametrisation of funSide closest to colPt
        s = fun.suppWidth/((Op.domain.component(funSide).dSv)*(colPt_R2-a_R2).');
        %compute distance to colPt from funSide
        d = norm(Op.domain.component(funSide).trace(s) - colPt_R2);
        %I = integral(F, a, b, 'arrayValued', true,'RelTol',1e-6);
        if s<=a
            [x,w] = grad_osc_quad(a,b,fun.suppWidth,d,Op.kwave,qppw,'L');
        elseif a<s && s<b
            [x1,w1] = grad_osc_quad(a,s,abs(s-a),d,Op.kwave,qppw,'R');
            [x2,w2] = grad_osc_quad(s,b,fun.suppWidth-abs(s-a),d,Op.kwave,qppw,'L');
            x = [x1; x2];
            w = [w1; w2];
        elseif b<=s
            [x,w] = grad_osc_quad(a,b,fun.suppWidth,d,Op.kwave,15,'R');
        else
            error('what is s then?');
        end
        I = w.'*F(x);
    end
end