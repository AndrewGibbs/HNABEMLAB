function I = bruteForce(Op, fun, funSide, colPt)
% use matlab's standard quadrature routine for integrals we haven't tackled
% appropriately yet:
    supp = fun.getSupp(funSide);
    a = supp(1);
    b = supp(2);
    F = @(y) Op.kernel(colPt.x, y, colPt.side, funSide) .* fun.eval(y);
    I = integral(F, a, b, 'arrayValued', true,'RelTol',1e-8);
end
