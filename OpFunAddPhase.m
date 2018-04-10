function G = OpFunAddPhase(Op,Fun,x, xGEy)
%for an operator S and function f with known phases, computes the phase of
%the integrand of Sf(x)

% ** hard code this for now:
OperatorDerivs = 3;
    
    %add common derivatives
    for n=1:min([OperatorDerivs+1 length(Fun.phase)])
        G{n} = @(y) Op.phaseAnalDeriv(x,y,n-1,xGEy) + Fun.phase{n}(y);%phaseAnalDeriv(s,t,n-1,xGEy)
    end
    
    %and fill in the rest
    if OperatorDerivs+1>length(Fun.phase)
        for n=length(Fun.phase)+1 : OperatorDerivs
            G{n} = @(y) Op.phaseAnalDeriv(x,y,n-1,xGEy);
        end
    elseif OperatorDerivs+1<length(Fun.phase)
        for n=OperatorDerivs+1 : length(Fun.phase)
            G{n} = @(y) Fun.phase{n}(y);
        end
    end
end

