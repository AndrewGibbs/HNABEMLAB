function G = OpFunAddPhase(Op,Fun,x, xGEy, NSDerivs)
%for an operator S and function f with known phases, computes the phase of
%the integrand of Sf(x)
    
    %add common derivatives
    for n=1:(NSDerivs+1)
        G{n} = @(y) Op.phaseAnalDeriv(x,y,n-1,xGEy) + Fun.phaseAnal(y,n-1,xGEy);%phaseAnalDeriv(s,t,n-1,xGEy)
    end
    
%     %and fill in the rest
%     if OperatorDerivs+1>length(Fun.phase)
%         for n=length(Fun.phase)+1 : OperatorDerivs
%             G{n} = @(y) Op.phaseAnalDeriv(x,y,n-1,xGEy);
%         end
%     elseif OperatorDerivs+1<length(Fun.phase)
%         for n=OperatorDerivs+1 : length(Fun.phase)
%             G{n} = @(y) Fun.phaseAnal(y,n-1,xGEy);
%         end
%     end
end

