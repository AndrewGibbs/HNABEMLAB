function [v, psi, ColMatrix, ColRHS] = ColHNA(Operator, HNAbasis, uinc, Gamma, varagin)
%computes oversampled collocation projection using HNA basis/frame
    Nquad = 15;
    %with RHS data
    f=DirichletData(uinc,Gamma);
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);
    
    DOFs=length(HNAbasis.el);
    
    %get collocation points:
    
    
    
    %** should eventually find a way to partition the basis into mesh
    %elements with + or - phase, so that quadrature can be reused.
    
    
    for n=1:length(X)
        logSingInfo=Singularity1D(X(n), Operator.singularity);
        for m=1:length(HNAbasis)
            
            %analytic extension of non-osc component of kernel:
            amp1 = @(y) Operator.kernelNonOscAnal(X(n),y) .* evalNonOscAnal(y);
            %and the corresponding phase:
            phase1 = OpFunAddPhase(Operator,HNAbasis.el(m),X(n));
            %now get weights and nodes:
            [ z1, w1 ] = NSD45( HNAbasis.el(m).a, HNAbasis.el(m).b, freq, Nquad, phase1,...
                        'fSingularities', logSingInfo, 'stationary points', [], 'order', []);
            %and evaluate integral:
            colMatrix(n,m) = (w1.'*amp1(z1));
        end
        phase2 = {@(y) sgn*(y-X(n)) + GOA.phaseLinear(1)*y + GOA.phaseLinear(2), @(y) sgn+GOA.phaseLinear(1), @(y) 0};
        [ z2, w2 ] = NSD45( HNAbasis.el(m).supp(1), HNAbasis.el(m).supp(2), freq, Nquad, phase2, 'fSingularities',logSingInfo);
        [~, valNonOsc]=GOA.eval(self,z2);
        %and evaluate the standard data and subtract the leading order integral:
        RHS(n) = f.eval(X) - w2.'*(amp2(z2).*valNonOsc);
    end
    
    %use least squares with Matlab's built in SVD to get coefficients
    coeffs=ColMatrix\ColRHS;
    v_N=Projection(coeffs,HNAbasis);
    
end