function [I, I2] = colEval(Op,fun, funSide, sCol, colSide, Nquad)

%function which evalutes integral Sf(x), essentially just a wrapper for
%NSD45 which ensures that the phase is always the analytic continuation of
%|x(s)-y(t)|

    %if function is defined over multiple sides, loop over these and sum up
    %contribution
    if length(funSide)>1
        I=0;
       for m=funSide
          I = I + colEval(Op,fun, m, sCol, colSide, Nquad); 
       end
       I2 = I;
       return;
    end

    %main function:
        maxSPorder = max(Op.phaseMaxStationaryPointOrder(funSide == colSide), fun.phaseMaxStationaryPointOrder);
        
        kwave = Op.kwave;
        
        %get endpoints of support of function
        supp = fun.getSupp(funSide);
        a = supp(1);
        b = supp(2);
        
        %return an error if we are this close to a singularity/branch point
        dangerZoneRad = 0.15*(b-a);
        singularSplit = dangerZoneRad;
    
    if funSide == colSide
        
        %same side singularity:
        distFun = @(t) abs(sCol - t);
        logSingInfo=singularity(sCol, Op.singularity, distFun);
        %choose the rectangle sufficiently small that phase is analytic
        rectrad = .5*min(logSingInfo.distFun(a),logSingInfo.distFun(b));
        
        %analytic extension of non-osc components of kernel:
        amp_a = @(y) Op.kernelNonOscAnal(sCol, y, true, colSide, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_b = @(y) Op.kernelNonOscAnal(sCol, y, false, colSide, funSide) .* fun.evalNonOscAnal(y, funSide);
        %and the corresponding phases:
        phase_a = OpFunAddPhase(Op, fun, funSide, sCol, colSide, true, maxSPorder+1);
        phase_b = OpFunAddPhase(Op, fun, funSide, sCol, colSide, false, maxSPorder+1);
        
        if  a < sCol && sCol < b
            %need to split the integral, as integrand not analytic at z=x

            if maxSPorder ==0
                [ z1a, w1a ] = NSD45( a, sCol, kwave, Nquad, phase_a,'settlerad',rectrad,...
                            'fSingularities', logSingInfo, 'stationary points', [], 'order', []);

                [ z1b, w1b ] = NSD45(sCol, b, kwave, Nquad, phase_b,'settlerad',rectrad,...
                            'fSingularities', logSingInfo, 'stationary points', [], 'order', []);
                I = (w1a.'*amp_a(z1a)) + (w1b.'*amp_b(z1b));

            else
                [ z1a, w1a ] = NSD45( a, sCol, kwave, Nquad, phase_a,'fSingularities', logSingInfo, 'settlerad', rectrad);

                [ z1b, w1b ] = NSD45(sCol, b, kwave, Nquad, phase_b,'fSingularities', logSingInfo, 'settlerad', rectrad);
                I = (w1a.'*amp_a(z1a)) + (w1b.'*amp_b(z1b));
            end
        else
            if sCol <= a
                amp = amp_b;
                phase = phase_b;
            elseif b <= sCol
                %analytic extension of non-osc component of kernel:
                amp = amp_a;
                phase = phase_a;
            else
                %this error will probably never ever happen:
                error('cant decide which is bigger of s and t');
            end
            %now get weights and nodes:
            if maxSPorder ==0
                [ z1, w1 ] = NSD45( a, b, kwave, Nquad, phase,...
                            'fSingularities', logSingInfo, 'stationary points', [], 'order', [], 'settlerad', rectrad);
            else
                [ z1, w1 ] = NSD45( a, b, kwave, Nquad, phase,'fSingularities', logSingInfo, 'settlerad', rectrad);
            end
            %and evaluate integral:
            I = (w1.'*amp(z1));
        end
        I2 = I;
    else %no branch in phase
        
        %different side singularity:
        distFun = @(t) Op.domain.distAnal(sCol, t, 0, [], colSide, funSide);
        %distR = Op.domain.distAnal(sCol, b, 0,[], colSide, funSide);
        logSingInfo=singularity([], Op.singularity, distFun);
        
        amp = @(y) Op.kernelNonOscAnal(sCol, y, [], colSide, funSide) .* fun.evalNonOscAnal(y, funSide);
        phase = OpFunAddPhase(Op, fun, funSide, sCol, colSide, [], maxSPorder+1);
        [stationaryPoints, orders, branchPoints] = symbolicStationaryPoints(Op.domain.side{colSide}.trace(sCol), fun, funSide, phase);
        [dangerTest, minCombo] =  min([abs(a-branchPoints(1)),abs(b-branchPoints(1)),abs(a-branchPoints(2)),abs(b-branchPoints(2))]);
            singularDifference = 0;
        if dangerTest < dangerZoneRad
            singularDifference = singularSplit;
            %warning('Highway to the dangerzone.');
            if ismember(minCombo,[1 3]) %singularity close to a
                a0 = a;
                a = a0 + singularSplit;
                [t, w0] = NonOsc45(a0, a, kwave, Nquad, phase{1}, logSingInfo, singularSplit);
                I0 = (w0.'*amp(t));
            elseif  ismember(minCombo,[2 4]) %singularity close to b
                b0 = b;
                b = b0 - singularSplit;
                [t, w0] = NonOsc45(b, b0, kwave, Nquad, phase{1}, logSingInfo, singularSplit);
                I0 = (w0.'*amp(t));
            end
        else
            I0 = 0;
        end
        rectrad = .5*min(logSingInfo.distFun(a),logSingInfo.distFun(b));
        %bodge this:
        if isa(fun,'GeometricalOpticsFunction')
           L =  fun.suppWidth(funSide) - singularDifference;
        else
            L = fun.L - singularDifference;
        end
        if isnan(stationaryPoints)
            %choose the rectangle sufficiently small that phase is analytic
            [ z, w ] = NSD45( a, b, kwave, Nquad, phase,'fSingularities', logSingInfo, 'settlerad', rectrad);
%             [X, W] = NonOsc45(a,b,kwave,Nquad,phase{1},logSingInfo, L);
        else %stationary points are already known
            [ z, w ] = NSD45( a, b, kwave, Nquad, phase,'fSingularities', logSingInfo, 'stationary points', stationaryPoints, 'order', orders, 'settlerad', rectrad);
%             [X, W] = NonOsc45(a,b,kwave,Nquad,phase{1},logSingInfo, L);
        end
        
        I = (w.'*amp(z)) + I0;
%         [X, W] = NonOsc45(a,b,kwave,Nquad,phase{1},logSingInfo, L);
%         I2 = (W.'*amp(X)) + I0;
        I2 = I;
        
        if isnan(I)
            warning('PathFinder returned a NaN, so using standard quadrature insteasd :-(');
            [X, W] = NonOsc45(a,b,kwave,Nquad,phase{1},logSingInfo, L);
            I = (W.'*amp(X)) + I0;
        end
        
    end
            
end