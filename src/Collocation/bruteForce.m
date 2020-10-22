function I = bruteForce(Op, fun, funSide, colPt)
% use matlab's standard quadrature routine for integrals we haven't tackled
% appropriately yet:
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
        I = integral(F, a, b, 'arrayValued', true,'RelTol',1e-8);
    end
end
