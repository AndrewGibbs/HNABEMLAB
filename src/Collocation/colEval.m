function [I, quadDataOut] = colEval(Op,fun, funSide, sCol, colSide, dist2a, dist2b, Nquad, quadDataIn, CGflag)

%function which evalutes integral Sf(x), essentially just a wrapper for
%NSD45 which ensures that the phase is always the analytic continuation of
%|x(s)-y(t)|

    roundErrorThresh = 1E-6;
    
    if nargin <= 9
        CGflag = false;
    end
    
    if CGflag
        minOscs = inf;
    else
        minOscs = 2;
    end
    %if function is defined over multiple sides, loop over these and sum up
    %contribution
    
    if length(funSide)>1
        I=0;
       for m=funSide
          I = I + colEval(Op,fun, m, sCol, dist2b, colSide, Nquad, CGflag); 
       end
       return;
    end
    if min(dist2a,dist2b)/abs(sCol)<roundErrorThresh
        b_roundErrorFlag = true;
    else
        b_roundErrorFlag = false;
    end
    
    %intiialise variables for data structure:
    z=[]; w=[]; z1a=[]; w1a=[]; z1b=[]; w1b=[];
     split = [0 0];

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
      
        %analytic extension of non-osc components of kernel:
        amp_a = @(y) Op.kernelNonOscAnal(sCol, y, true, colSide, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_b = @(y) Op.kernelNonOscAnal(sCol, y, false, colSide, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_a_flip = @(r) Op.kernelNonOscAnal(r, 0, true, colSide, funSide).* fun.evalNonOscAnal(sCol - r, funSide);
        amp_b_flip = @(r) Op.kernelNonOscAnal(r, 0, true, colSide, funSide).* fun.evalNonOscAnal(sCol + r, funSide);
        %and the corresponding phases:
        phase_a = OpFunAddPhase(Op, fun, funSide, sCol, colSide, true, maxSPorder+1);
        phase_b = OpFunAddPhase(Op, fun, funSide, sCol, colSide, false, maxSPorder+1);
        %phase_b_flip = OpFunAddPhase(Op, fun, funSide, sCol, colSide, false, maxSPorder+1);
        
        for n = 1:length(phase_b)
            phase_a_flip{n} = @(r) (-1)^(n+1)*phase_a{n}(sCol-r);
            phase_b_flip{n} = @(r) phase_b{n}(r+sCol);
        end
        
        %now the more general amp, for when there is no branch in [a,b]
        amp = @(y) Op.kernelNonOscAnal(sCol, y, [], colSide, funSide) .* fun.evalNonOscAnal(y, funSide);
        phase = OpFunAddPhase(Op, fun, funSide, sCol, colSide, [], maxSPorder+1);
        
    if ~isempty(quadDataIn)
        if isequal(quadDataIn.split,[0 0])
            I = (quadDataIn.w.'*amp(quadDataIn.z));
        else
            I = 0;
            %add components seperately
            if quadDataIn.split(1) == 1
                I = I + quadDataIn.w1a.'*amp_a(quadDataIn.z1a);
            end
            if quadDataIn.split(2) == 1
                I = I + quadDataIn.w1b.'*amp_b(quadDataIn.z1b);
            end
        end
        quadDataOut = quadDataIn;
        return;
    end
        
    if funSide == colSide
        
        %same side singularity:
        distFun = @(t) abs(sCol - t);
        logSingInfo=singularity(sCol, Op.singularity, distFun);
        %choose the rectangle sufficiently small that phase is analytic
        rectrad = .5*min(logSingInfo.distFun(a),logSingInfo.distFun(b));
        
        if  a < sCol && sCol < b
            %need to split the integral, as integrand not analytic at z=x

            if maxSPorder ==0

                if b_roundErrorFlag
                    
                    logSingInfo_flip_a = logSingInfo;
                    logSingInfo_flip_a.position = 0;
                    logSingInfo_flip_a.distFun = @(r) abs(r);
                    [ z1a, w1a ] = PathFinder( 0, dist2a, kwave, Nquad, phase_a_flip,'settlerad',rectrad,...
                                'fSingularities', logSingInfo_flip_a, 'stationary points', [], 'order', [],'minOscs',minOscs);
                     I1 = (w1a.'*amp_a_flip(z1a));
                 
                    logSingInfo_flip_b = logSingInfo;
                    logSingInfo_flip_b.position = 0;
                    logSingInfo_flip_b.distFun = @(r) abs(r);
                    [ z1b, w1b ] = PathFinder(0, dist2b, kwave, Nquad, phase_b_flip,'settlerad',rectrad,...
                                'fSingularities', logSingInfo_flip_b, 'stationary points', [], 'order', [],'minOscs',minOscs);
                    I2 = (w1b.'*amp_b_flip(z1b));
                    1+1;
                else
                    
                    [ z1a, w1a ] = PathFinder( a, sCol, kwave, Nquad, phase_a,'settlerad',rectrad,...
                                'fSingularities', logSingInfo, 'stationary points', [], 'order', [],'minOscs',minOscs);
                     I1 = (w1a.'*amp_a(z1a));
                 
                    [ z1b, w1b ] = PathFinder(sCol, b, kwave, Nquad, phase_b,'settlerad',rectrad,...
                                'fSingularities', logSingInfo, 'stationary points', [], 'order', [],'minOscs',minOscs);
                    I2 = (w1b.'*amp_b(z1b));
                end
                
                I = I1 + I2;

            else
                [ z1a, w1a ] = PathFinder( a, sCol, kwave, Nquad, phase_a,'fSingularities', logSingInfo, 'settlerad', rectrad,'minOscs',minOscs);
                I1 = (w1a.'*amp_a(z1a)) ;
                
                [ z1b, w1b ] = PathFinder(sCol, b, kwave, Nquad, phase_b,'fSingularities', logSingInfo, 'settlerad', rectrad,'minOscs',minOscs);
                I2 = (w1b.'*amp_b(z1b));

%                 if b_roundErrorFlag
%                     [ z1b, w1b ] = PathFinder(0, dist2b, kwave, Nquad, phase_b_flip,'fSingularities', logSingInfo, 'settlerad', rectrad,'minOscs',minOscs);
%                     I2 = (w1b.'*amp_b_flip(z1b));
%                 else
%                     [ z1b, w1b ] = PathFinder(sCol, b, kwave, Nquad, phase_b,'fSingularities', logSingInfo, 'settlerad', rectrad,'minOscs',minOscs);
%                     I2 = (w1b.'*amp_b(z1b));
%                 end
                
                I = I1 + I2;
            end
            
            split = [1 1];
        else
            if sCol <= a
                amp = amp_b;
                phase = phase_b;
                phase_flip = phase_b_flip;
                split = [0 1];
                dist2end = dist2b;
            elseif b <= sCol
                %analytic extension of non-osc component of kernel:
                amp = amp_a;
                phase = phase_a;
                phase_flip = phase_a_flip;
                split = [1 0];
                dist2end = dist2a;
            else
                %this error will probably never ever happen:
                error('cant decide which is bigger of s and t');
            end
            %now get weights and nodes:
            if maxSPorder ==0
                [ z_, w_ ] = PathFinder( a, b, kwave, Nquad, phase,...
                            'fSingularities', logSingInfo, 'stationary points', [], 'order', [], 'settlerad', rectrad,'minOscs',minOscs,'width',fun.suppWidth);
                        
%                 logSingInfo_flip = logSingInfo;
%                 logSingInfo_flip.position = 0;
%                 logSingInfo_flip.distFun = @(r) abs(r);
%                 [ z_, w_ ] = PathFinder( 0, fun.suppWidth, kwave, Nquad, phase_flip,'settlerad',rectrad,...
%                         'fSingularities', logSingInfo_flip, 'stationary points', [], 'order', [], 'minOscs', minOscs);
            else
                [ z_, w_ ] = PathFinder( a, b, kwave, Nquad, phase,'fSingularities', logSingInfo, 'settlerad', rectrad,'minOscs',minOscs,'width',fun.suppWidth);
            end
            %and evaluate integral:
            I = w_.'*amp(z_);
            %now store in correct form:
            if sCol <= a
                w1b = w_;
                z1b = z_;
            else
                w1a = w_;
                z1a = z_;
            end
        end
    else %no branch in phase
        
        %different side singularity:
        distFun = @(t) Op.domain.distAnal(sCol, t, 0, [], colSide, funSide);
        %distR = Op.domain.distAnal(sCol, b, 0,[], colSide, funSide);
        logSingInfo=singularity([], Op.singularity, distFun);
        
        [stationaryPoints, orders, branchPoints] = symbolicStationaryPoints(Op.domain.side{colSide}.trace(sCol), fun, funSide, phase);
        [dangerTest, minCombo] =  min([abs(a-branchPoints(1)),abs(b-branchPoints(1)),abs(a-branchPoints(2)),abs(b-branchPoints(2))]);
            singularDifference = 0;
        if dangerTest < dangerZoneRad
            singularDifference = singularSplit;
            if ismember(minCombo,[1 3]) %singularity close to a
                a0 = a;
                a = a0 + singularSplit;
                [t, w0] = NonOsc45(a0, a, kwave, Nquad, phase{1}, logSingInfo, singularSplit);
                1+1;
                %I0 = (w0.'*amp(t));
            elseif  ismember(minCombo,[2 4]) %singularity close to b
                b0 = b;
                b = b0 - singularSplit;
                [t, w0] = NonOsc45(b, b0, kwave, Nquad, phase{1}, logSingInfo, singularSplit);
                1+1;
                %I0 = (w0.'*amp(t));
            end
        else
            w0 = []; t = [];
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
            [ z1, w1 ] = PathFinder( a, b, kwave, Nquad, phase,'fSingularities', logSingInfo, 'settlerad', rectrad,'minOscs',minOscs);
        else %stationary points are already known
            [ z1, w1 ] = PathFinder( a, b, kwave, Nquad, phase,'fSingularities', logSingInfo, 'stationary points', stationaryPoints, 'order', orders, 'settlerad', rectrad,'minOscs',minOscs);
        end
        z = [t; z1];   w = [w0; w1];
        I = (w.'*amp(z));
        
        if isnan(I)
            warning('PathFinder returned a NaN, so using standard quadrature insteasd :-(');
            [X, W] = NonOsc45(a,b,kwave,Nquad,phase{1},logSingInfo, L);
            z = [X; t];
            W = [W; w0];
            I = (W.'*amp(X)) + (w0.'*amp(t));
        end
        
    end
    
    quadDataOut = struct('z',z,'w',w,'split',split,'z1a', z1a, 'w1a',w1a,'z1b', z1b, 'w1b', w1b);
            
end