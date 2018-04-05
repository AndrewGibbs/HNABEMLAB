function G = OpFunAddPhase(Op,Fun,x)
%for an operator S and function f with known phases, computes the phase of
%the integrand of Sf(x)
    
    %add common derivatives
    for n=min(length(Op.phase),length(Fun.phase))
        G{n} = @(y) Op.phase{n}(x,y) + Fun.phase{n}(y);
    end
    
    %and fill in the rest
    if length(Op.phase)>length(Fun.phase)
        for n=length(Fun.phase)+1 : length(Op.phase)
            G{n} = @(y) Op.phase{n}(x,y);
        end
    elseif length(Op.phase)<length(Fun.phase)
        for n=length(Op.phase)+1 : length(Fun.phase)
            G{n} = @(y) Fun.phase{n}(y);
        end
    end
end

