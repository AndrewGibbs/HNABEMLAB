function G = OpFunAddPhaseCorner(Op, Fun, funSide, colDist, supp2corner, colSide, NSDerivs ,D, flip)
%for an operator S and function f with known phases, computes the phase of
%the integrand of Sf(x)
    
    %add common derivatives
    for n=1:(NSDerivs+1)
        if ~flip
            G{n} = @(t) Op.phaseAnalDerivCorner(colDist, t, supp2corner, n-1, colSide, funSide)...
                   + Fun.phaseAnal( D+ t, n-1, funSide);%+ Fun.phaseAnal( supp2corner + t, n-1, funSide);
        else
            G{n} = @(t) Op.phaseAnalDerivCorner(colDist, t, supp2corner, n-1, colSide, funSide)...
                   + (-1)^(n-1)*Fun.phaseAnal( D - t, n-1, funSide); %*Fun.phaseAnal( supp2corner - t, n-1, funSide);
        end
    end
    
end