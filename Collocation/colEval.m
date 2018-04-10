function I = colEval(Op,fun,x)

%function which evalutes integral Sf(x), essentially just a wrapper for
%NSD45 which ensures that the phase is always the analytic continuation of
%|x(s)-y(t)|
%S is an integral operator here

        logSingInfo=Singularity1D(x, Op.singularity);
        kwave = Op.kwave;
        Nquad = 15;
        
        %analytic extension of non-osc components of kernel:
        amp_a = @(y) Op.kernelNonOscAnal(x,y, true) .* fun.evalNonOscAnal(y);
        amp_b = @(y) Op.kernelNonOscAnal(x,y, false) .* fun.evalNonOscAnal(y);
        %and the corresponding phases:
        phase_a = OpFunAddPhase(Op,fun,x,true);
        phase_b = OpFunAddPhase(Op,fun,x,false);
        
    if  fun.a < x && x < fun.b
        %need to split the integral, as integrand not analytic at z=x

        [ z1a, w1a ] = NSD45( fun.a, x, kwave, Nquad, phase_a,...
                    'fSingularities', logSingInfo, 'stationary points', [], 'order', []);

        [ z1b, w1b ] = NSD45(x, fun.b, kwave, Nquad, phase_b,...
                    'fSingularities', logSingInfo, 'stationary points', [], 'order', []);
        I = (w1a.'*amp_a(z1a)) + (w1b.'*amp_b(z1b));
    else
        if x <= fun.a
            amp = amp_b;
            phase = phase_b;
        elseif fun.b <= x
            %analytic extension of non-osc component of kernel:
            amp = amp_a;
            phase = phase_a;
        else
            %this error will probably never ever happen:
            error('cant decide which is bigger of s and t');
        end
        %now get weights and nodes:
        [ z1, w1 ] = NSD45( fun.a, fun.b, kwave, Nquad, phase,...
                    'fSingularities', logSingInfo, 'stationary points', [], 'order', []);
        %and evaluate integral:
        I = (w1.'*amp(z1));
    end
            
            
end

