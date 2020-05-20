function [I, quadDataOut] = LHSquad(Op,fun, funSide, colPt, Nquad, quadDataIn)
%decides which quadrature rule to use for matrix collocation entries
    if colPt.side == funSide
        % same as screen, can use singular NSD here
        [I, quadDataOut] = colEvalLinear(Op,fun, funSide, colPt, Nquad, quadDataIn);
    else
        % more complex, use brute force quadrature for now
        I = bruteForce(Op,fun, funSide, colPt);
        quadDataOut = [];
    end
end

