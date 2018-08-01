function G = OpFunAddPhase(Op, Fun, funSide, s, sSide, sGEt, NSDerivs)
%for an operator S and function f with known phases, computes the phase of
%the integrand of Sf(x)
    
    %add common derivatives
    for n=1:(NSDerivs+1)
        G{n} = @(t) Op.phaseAnalDeriv(s,t,n-1,sGEt,sSide,funSide) + Fun.phaseAnal(t,n-1,funSide);
    end
    
end