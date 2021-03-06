function [I, quadDataOut] = colEvalV4(Op,fun, funSide, colPt, Nquad, quadDataIn, CGflag, linearSpecialCase)
%linearSpecialCase = false;
%function which evalutes integral Sf(x), essentially just a wrapper for
%NSD45 which ensures that the phase is always the analytic continuation of
%|x(s)-y(t)|

%sCol, colSide, dist2a, dist2b - all absorbed into X
    
    if nargin <= 6
        CGflag = false;
    end
    
    if CGflag
        minOscs = inf;
    else
        minOscs = 5;
    end
    %if function is defined over multiple sides, loop over these and sum up
    %contribution
    
    if length(funSide)>1
        I=0;
       for m=funSide
           if ~isempty(quadDataIn)
                [I_, ~] = colEvalV2(Op,fun, m, colPt.x, dist2b, colSide, Nquad, quadDataIn{m}, CGflag); 
           else
               [I_, quadDataOut{m}] = colEvalV2(Op,fun, m, colPt.x, dist2b, colSide, Nquad, [], CGflag); 
           end
          I = I + I_;
       end
       return;
    end
    
%     %intiialise variables for data structure:
%     z=[]; w=[]; z1a=[]; w1a=[]; z1b=[]; w1b=[];
%      split = [0 0];

    %main function:
        maxSPorder = max(Op.phaseMaxStationaryPointOrder(funSide == colPt.side), fun.phaseMaxStationaryPointOrder);
        
        kwave = Op.kwave;
        
        %get endpoints of support of function
        supp = fun.getSupp(funSide);
        a = supp(1);
        b = supp(2);
        
        %return an error if we are this close to a singularity/branch point
        dangerZoneRad = .35;%0.25/kwave;%max(0.15*(b-a),dangerWidth);
        %singularSplit = dangerZoneRad;
        
        dangerWidth = 1E-12;
      
        %analytic extension of non-osc components of kernel:
        amp_a = @(y) Op.kernelNonOscAnal(colPt.x, y, true, colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_b = @(y) Op.kernelNonOscAnal(colPt.x, y, false, colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_a_flip = @(r) Op.kernelNonOscAnal(r, 0, true, colPt.side, funSide).* fun.evalNonOscAnal(colPt.x - r, funSide);
        amp_b_flip = @(r) Op.kernelNonOscAnal(r, 0, true, colPt.side, funSide).* fun.evalNonOscAnal(colPt.x + r, funSide);
        %and the corresponding phases:
        phase_a = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, true, maxSPorder+1);
        phase_b = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, false, maxSPorder+1);
        %phase_b_flip = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, false, maxSPorder+1);
        
        for n = 1:length(phase_b)
            phase_a_flip{n} = @(r) (-1)^(n+1)*phase_a{n}(colPt.x-r);
            phase_b_flip{n} = @(r) phase_b{n}(r+colPt.x);
        end
        
        %now the more general amp, for when there is no branch in [a,b]
        amp = @(y) Op.kernelNonOscAnal(colPt.x, y, [], colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide);
        phase = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, [], maxSPorder+1);
        
        
        if length(fun.domain.L )>1
            L = fun.domain.L(funSide);
        else
            L = fun.domain.L;
        end
            
    if funSide == colPt.side
        %same side singularity:
        distFun = @(t) abs(colPt.x - t);
        logSingInfo=singularity(colPt.x, Op.singularity, distFun);
        %choose the rectangle sufficiently small that phase is analytic
        rectrad = .5*min(logSingInfo.distFun(a),logSingInfo.distFun(b));
        
        if maxSPorder ==0
            SPin = [];
            SPOin = [];
        else
            SPin = NaN;
            SPOin = NaN;
        end
        
        if ((fun.meshEl.distR < colPt.distSideR) && (fun.meshEl.distR + fun.meshEl.width > colPt.distSideR) && (colPt.distSideR <= colPt.distSideL)) ...
             || ((fun.meshEl.distL < colPt.distSideL) && (fun.meshEl.distL + fun.meshEl.width > colPt.distSideL) && (colPt.distSideL <= colPt.distSideR))
            %the above is a rounding-error-friendly version of:
            % a < colPt.x && colPt.x < b
           
            %need to split the integral, as integrand not analytic at z=x
            
            if isempty(colPt.distMeshL) || isempty(colPt.distMeshR)  %need to bodge this, for the overlapping case only
                %cannot determine which mesh element we're on, given
                %collocation point
                if colPt.distSideL < colPt.distSideR
                    colPt.distMeshL = colPt.x - fun.meshEl.interval(1);
                    colPt.distMeshR = fun.meshEl.interval(2) - colPt.x;
                else
                    colPt.distMeshL = colPt.distSideL - fun.meshEl.distL;
                    colPt.distMeshR = colPt.distSideR - fun.meshEl.distR;
                end
                    %colPt.distMeshR = colPt.distSideR - fun.meshEl.distR;
                %end
            end

           if ~isempty(quadDataIn)
                w1a = quadDataIn.w1a;
                w1b = quadDataIn.w1b;
                z1a = quadDataIn.z1a;
                z1b = quadDataIn.z1b;
            else
                logSingInfo_flip_a = logSingInfo;
                logSingInfo_flip_a.position = 0;
                logSingInfo_flip_a.distFun = @(r) abs(r);
                if linearSpecialCase
                    alphaPhase_a = phase_a_flip{2}(1);
                    missingPhaseShift_a = phase_a_flip{1}(1) - alphaPhase_a;
                    [ z1a, w1a ] = LinearNSD(colPt.distMeshL, kwave, alphaPhase_a, 0, Nquad);
                    w1a = w1a*exp(1i*kwave*missingPhaseShift_a);
                else
                    [ z1a, w1a ] = PathFinder( 0, colPt.distMeshL, kwave, Nquad, phase_a_flip,'settlerad',rectrad,...
                            'fSingularities', logSingInfo_flip_a, 'stationary points', SPin, 'order', SPOin,'minOscs',minOscs,'linear');
                end

                logSingInfo_flip_b = logSingInfo;
                logSingInfo_flip_b.position = 0;
                logSingInfo_flip_b.distFun = @(r) abs(r);
                
                if linearSpecialCase
                    alphaPhase_b = phase_b_flip{2}(1);
                    missingPhaseShift_b = phase_b_flip{1}(1) - alphaPhase_b;
                    [ z1b, w1b ] = LinearNSD(colPt.distMeshR, kwave, alphaPhase_b, 0, Nquad);
                    w1b = w1b*exp(1i*kwave*missingPhaseShift_b);
                else
                    [ z1b, w1b ] = PathFinder(0, colPt.distMeshR, kwave, Nquad, phase_b_flip,'settlerad',rectrad,...
                            'fSingularities', logSingInfo_flip_b, 'stationary points', SPin, 'order', SPOin,'minOscs',minOscs,'linear');
                end
                    
                quadDataOut.w1a = w1a;
                quadDataOut.w1b = w1b;
                quadDataOut.z1a = z1a;
                quadDataOut.z1b = z1b;
           end
                if sum(abs(w1a))>0
                    I1 = (w1a.'*amp_a_flip(z1a));
                else
                    I1 = 0;
                end
                if sum(abs(w1b))>0
                    I2 = (w1b.'*amp_b_flip(z1b));
                else
                    I2 = 0;
                end
                I = I1 + I2;
        else %same side, no singularity in (a,b)
            
            if colPt.distSideL>colPt.distSideR %col pt is closer to right endpoint of mesh element than left
                %so deal in terms of distance from right endpoint
                if colPt.distSideR<=fun.meshEl.distR % a<b<x
                    type_ab = 'b';
                    phase = phase_a;
                else % x<a<b
                    type_ab = 'a';
                    phase = phase_b;
                end
            else %otherwise deal in terms of distance from left endpoint
                if colPt.distSideL<=fun.meshEl.distL % x<a<b
                    type_ab = 'a';
                    phase = phase_b;
                else % a<b<x
                    type_ab = 'b';
                    phase = phase_a;
                end
            end
            
%             if colPt.x <= a
%                 type_ab = 'a';
%                 phase = phase_b;
%             elseif b <= colPt.x
%                 type_ab = 'b';
%                 phase = phase_a;
%             else
%                 %this error will probably never ever happen:
%                 error('cant decide which is bigger of s and t');
%             end
            
            if fun.meshEl.distL <= fun.meshEl.distR
                type_LR = 'L';
                a_star = a;
                b_star = b;
                sing_star_point = colPt.x;
                phase_star = phase;
            else
                type_LR = 'R';
                a_star = fun.meshEl.distR; %= L-b;
                b_star = fun.meshEl.distR + fun.meshEl.width;  %= L-a;
                L_minus_a = b_star;
                sing_star_point = colPt.distSideR;
                for n = 1:length(phase_b)
                    phase_star{n} = @(z) (-1)^(n+1)*phase{n}(L-z);
                end
            end
            
            %create singularity info:
            logSingInfo_star = singularity(sing_star_point, Op.singularity);
            
            %there are four cases for the amp
            switch strcat(type_LR,type_ab)
                case 'La'
                    amp_star = amp_b;
                case 'Lb'
                    amp_star = amp_a;
                case 'Ra'
                    amp_star =  @(z) Op.kernelNonOscAnal(colPt.distSideR, z, true, colPt.side, funSide) .* fun.evalNonOscAnalPivot(z, funSide, L_minus_a);
                case 'Rb'
                    amp_star =  @(z) Op.kernelNonOscAnal(colPt.distSideR, z, false, colPt.side, funSide) .* fun.evalNonOscAnalPivot(z, funSide, L_minus_a);
            end
            
            if ~isempty(quadDataIn)
                w_ = quadDataIn.w_;
                z_ = quadDataIn.z_;
            else
                %now get weights and nodes:
                if linearSpecialCase
                    T = fun.meshEl.width;
                    if type_ab == 'b'
                        if type_LR == 'L'
                            xi = fun.supp(2) - colPt.x;
                        else
                            xi = colPt.distSideR-(fun.meshEl.distR);
                        end
                    else
                        if type_LR == 'L'
                            xi =  colPt.x-fun.supp(1);
                        else
                            xi = (fun.meshEl.distR + fun.meshEl.width) - colPt.distSideR ;%colPt.distSideR-fun.meshEl.distR;
                        end
                    end
                    %phaseSign = (strcmp(phase_star{2}(1),'b')*-2+1)*(strcmp(type_LR,'R')*-2+1);
                    %
                    
                    if phase_star{2}(1)<0
                        toBeFlipped = true;
                        alphaPhase = -phase_star{2}(1);
                    else
                        toBeFlipped = false;
                        alphaPhase = phase_star{2}(1);
                    end
                    
                    if type_ab == 'b'
                        if type_LR == 'L'
                            missingPhaseShift = phase_star{1}(1) - phase_star{2}(1) - fun.supp(2)*alphaPhase;
                        else
                            %missingPhaseShift = -(phase_star{1}(1) - phase_star{2}(1)) - fun.supp(2)*alphaPhase;%!
                            missingPhaseShift = phase_star{1}(0) + a_star*alphaPhase;%!
                        end
                    else
                        if type_LR == 'L'
                            missingPhaseShift = phase_star{1}(1) - phase_star{2}(1) + fun.supp(1)*alphaPhase;
                        else
                            missingPhaseShift = phase_star{1}(0) - b_star*alphaPhase;
                        end
                    end
                    
                    [ z_, w_ ] = LinearNSD(T, kwave, alphaPhase, xi, Nquad);
                    w_ = w_*exp(1i*kwave*missingPhaseShift);
                    if toBeFlipped
                       z_ = b_star-z_;
                    else
                       z_ = a_star+z_;
                    end
                    
%                     if type_ab == 'a'
%                         amp_star = @(y) Op.kernelNonOscAnal(0, y+abs(xi), true, colPt.side, funSide) .* fun.evalNonOscAnal(y+fun.supp(1), funSide);
%                     else
%                         amp_star = @(y) Op.kernelNonOscAnal(0, fun.supp(2)-y-colPt.x, true, colPt.side, funSide) .* fun.evalNonOscAnal(fun.supp(2)-y, funSide);
%                     end
                else
                    [ z_, w_ ] = PathFinder( a_star, b_star, kwave, Nquad, phase_star,'settlerad',rectrad,...
                            'fSingularities', logSingInfo_star, 'stationary points', SPin, 'order', SPOin, 'minOscs', minOscs, 'width', fun.suppWidth,'linear');
                end
                    
                quadDataOut.w_ = w_;
                quadDataOut.z_ = z_;
            end
                
            %and evaluate integral:
            I = w_.'*amp_star(z_);
            %now store in correct form:
            if colPt.x <= a
                w1b = w_;
                z1b = z_;
            else
                w1a = w_;
                z1a = z_;
            end
        end
    elseif isa(Op.domain,'MultiScreen')
        
        if colPt.side > funSide % x_1(s) > Y_1(t)
            shiftedCol = Op.domain.segSplits(2*colPt.side-1) - Op.domain.segSplits(2*funSide-1)  +  colPt.x;
            amp_shift = @(y) Op.kernelNonOscAnal(shiftedCol, y, true, colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide); 
            phase_shift = OpFunAddPhase(Op, fun, funSide, shiftedCol, colPt.side, true, maxSPorder+1);
        else                    % x_1(s) < Y_1(t)
            shiftedCol = -(Op.domain.segSplits(2*funSide-1) - Op.domain.segSplits(2*colPt.side-1))  +  colPt.x;
            amp_shift = @(y) Op.kernelNonOscAnal(shiftedCol, y, false, colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide); 
            phase_shift = OpFunAddPhase(Op, fun, funSide, shiftedCol, colPt.side, false, maxSPorder+1);
        end
        distFun = @(t) abs(shiftedCol - t);
        logSingInfo = singularity(shiftedCol, Op.singularity, distFun);
        rectRad = .5*min(logSingInfo.distFun(a),logSingInfo.distFun(b));
        
       if ~isempty(quadDataIn)
            w_ = quadDataIn.w_;
            z_ = quadDataIn.z_;
        else
            %now get weights and nodes:
            if a==b %no singularity & zero width, contribution will be negligable
                z_ = a;
                w_ = 0;
            else
                
                    
            if false %unfinished - need to fix rounding errors in (e.g.) "betaConst-alphaPhase*b",
                    % when these are all large, get different values to
                    % PathFinder and integral(.,.,.). However, the abs
                    % value of the weights is the same, suggesting that
                    % the error is only a rotation in the complex
                    % plane.
                alphaPhase = phase_shift{2}(0);
                betaConst = phase_shift{1}(0);
                missingPhaseShift = phase_shift{1}(1) - alphaPhase + a;

                if alphaPhase>0
                    [ z_, w_ ] = LinearNSD(fun.suppWidth, kwave, alphaPhase, shiftedCol, Nquad);
                    %w_ = w_*exp(1i*kwave*(betaConst+a*alphaPhase));
                    w_ = w_*exp(1i*kwave*betaConst)*exp(1i*kwave*a*alphaPhase);
                else
                    [ z_, w_ ] = LinearNSD(fun.suppWidth, kwave, -alphaPhase, -shiftedCol, Nquad);
                    w_ = w_*exp(1i*kwave*(betaConst-alphaPhase*b));
                    %w_ = w_* exp(1i*kwave*betaConst)*exp(-1i*kwave*alphaPhase*b);
                    z_ = b-z_;
                end
            else
                [ z_, w_ ] = PathFinder( a, b, kwave, Nquad, phase_shift,'settlerad', rectRad, 'fSingularities', logSingInfo,...
                                        'stationary points', [], 'order', [], 'minOscs', minOscs, 'width', fun.suppWidth,'linear');
            end
        end
        quadDataOut.w_ = w_;
        quadDataOut.z_ = z_;
        end

        %and evaluate integral:
        I = w_.'*amp_shift(z_);
            
%        I_check = integral(@(t) amp_shift(t).*exp(1i*kwave.*phase_shift{1}(t)),a,b,'ArrayValued',true);
%         quadDataOut = [];
        
    else %no branch in phase
        
%         I = integral(@(z) amp(z).*exp(1i*kwave*phase{1}(z)), a, b, 'arrayValued', true,'RelTol',1e-13);
%         quadDataOut = [];
%         return;

        quadDataOut = [];
        %[stationaryPoints, orders, branchPoints] = symbolicStationaryPoints(Op.domain.component(colPt.side).trace(colPt.x), fun, funSide, phase);
        stationaryPoints = [];
        orders = [];
        branchPoints = [];
        
        SPin = stationaryPoints;
        SPOin = orders;
        
        %determine if endpoint is close to singularity on neighbouring
        %side:
        if ~isempty(SPin)
            if min(abs(a-branchPoints(1)),abs(a-branchPoints(2))) < dangerZoneRad
                singClose2a = true;
            else
                singClose2a = false;
            end

            if min(abs(b-branchPoints(1)),abs(b-branchPoints(2))) < dangerZoneRad
                singClose2b = true;
            else
                singClose2b = false;
            end

            grad2a = false; grad2b = false;
            if singClose2a
                if a>L/2
                    grad2b = true;
                else
                    grad2a = true;
                end
            elseif singClose2b
                if b<L/2
                    grad2a = true;
                else
                    grad2b = true;
                end
            end
        else
            grad2a = false;
            grad2b = false;
        end
        
        if ~grad2a && ~ grad2b %no singularity issues
            
            if ~isempty(quadDataIn)
                %load quad data and skip to the sum
                z = quadDataIn.z;
                w = quadDataIn.w;
            else
                if isa(fun,'GeometricalOpticsEdge')
                    width = fun.suppWidth;
                else
                    width = fun.suppWidth;
                end
                distFun = @(t) Op.domain.distAnal(colPt.x, t, 0, [], colPt.side, funSide);
                %distR = Op.domain.distAnal(colPt.x, b, 0,[], colPt.side, funSide);
                logSingInfo=singularity([], Op.singularity, distFun);
                rectrad = .5*min(logSingInfo.distFun(a),logSingInfo.distFun(b));
                if width<minOscs*2*pi/kwave
                    [ z, w ] = PathFinder( a, b, kwave, Nquad, phase, ...
                                            'stationary points', stationaryPoints, 'order', orders, 'settlerad', ...
                                                rectrad,'minOscs',inf, 'width', width);
                                            %have changed minOscs to inf, to
                                            %always use standard quad here ^^
                    %I = (w.'*amp(z));
                else

                    if ~isempty(stationaryPoints)
                        if  ~(a <= stationaryPoints && stationaryPoints <= b)
                            stationaryPoints = [];
                            orders = [];
                        end
                    end
                    if isempty(stationaryPoints)
                        %perhaps the phase is flat near the endpoints, so chop
                        %off the first few oscillations:
                        XoscL = findNonOscBit(phase{1},a,b,kwave,minOscs);
                        XoscR = findNonOscBitR(phase{1},a,b,kwave,minOscs);
                        if XoscR<XoscL
                            %do the whole thing non-oscilllatorily
                            [ z, w ] = PathFinder( a, b, kwave, Nquad, phase, ...
                                            'stationary points', stationaryPoints, 'order', orders, 'settlerad', ...
                                                rectrad,'minOscs',inf, 'width', width);
                        else
                            [ za, wa ] = PathFinder( a, XoscL, kwave, Nquad, phase, ...
                                            'stationary points', stationaryPoints, 'order', orders, 'settlerad', ...
                                                rectrad,'minOscs',inf, 'width', width);

                            [ z_mid, w_mid ] = PathFinder( XoscL, XoscR, kwave, Nquad, phase, ...
                                            'stationary points', stationaryPoints, 'order', orders, 'settlerad', ...
                                                rectrad,'minOscs',minOscs, 'width', width);
                            [ zb, wb ] = PathFinder( XoscR, b, kwave, Nquad, phase, ...
                                            'stationary points', stationaryPoints, 'order', orders, 'settlerad', ...
                                               rectrad,'minOscs',inf, 'width', width);
                             z = [za; z_mid; zb];
                             w = [wa; w_mid; wb];
                        end
                        %I = PathFinderChebWrap(a,b,kwave,Nquad,amp,phase,logSingInfo,stationaryPoints);
                    else% a <= stationaryPoints && stationaryPoints <= b
                         [ z, w ] = PathFinder( a, b, kwave, Nquad, phase, ...
                                            'stationary points', stationaryPoints, 'order', orders, 'settlerad', ...
                                                rectrad,'minOscs',minOscs, 'width', width);
                    end
                end
                %save quad data
                quadDataOut.z = z;
                quadDataOut.w = w;
            end
            
            I = (w.'*amp(z));
            return;
        end
        
        if grad2a
            if isa(fun,'GeometricalOpticsFunction')
                a_shift = 0;
                b_shift = fun.suppWidth(funSide);
                width = fun.suppWidth(funSide);
            else
                a_shift = fun.meshEl.distL; %want dista
                b_shift = fun.meshEl.distL + fun.suppWidth;
                width = fun.suppWidth;
            end
            suppDistCorner = a_shift;
            
            colDistCorner = colPt.distSideR;
            
            amp_corner = @(y) Op.kernelNonOscAnalCorner( colPt.distSideR, y, a_shift, colPt.side, funSide) .* fun.evalNonOscAnal( a + y, funSide); %changed from (a_shift +y, funSide);
            phaseCorner = OpFunAddPhaseCorner(Op, fun, funSide, colPt.distSideR, a_shift, colPt.side, maxSPorder+1 , a, false);
            
            
        elseif grad2b
            if isa(fun,'GeometricalOpticsFunction')
                a_shift = 0;
                b_shift = fun.suppWidth(funSide);
                width = fun.suppWidth(funSide);
            else
                a_shift = fun.meshEl.distR;
                b_shift = fun.meshEl.distR + fun.suppWidth;
                width = fun.suppWidth;
            end
            suppDistCorner = a_shift;
            
            colDistCorner = colPt.distSideL;
            
            amp_corner = @(y) Op.kernelNonOscAnalCorner( colPt.distSideL, y, a_shift, colPt.side, funSide).* fun.evalNonOscAnal( b - y, funSide);  %changed from (a_shift -y, funSide);
            phaseCorner = OpFunAddPhaseCorner(Op, fun, funSide, colPt.distSideL, a_shift, colPt.side, maxSPorder+1 , b, true);
        end
        
            %determine the location of the singularities after this change
            %of variables
            internalAngle = Op.domain.internalAngle(colPt.side,funSide);
            %sing_flip = mean(roots([1,2*suppDistCorner - 2*colDistCorner*cos(internalAngle),suppDistCorner^2 + colDistCorner^2 -2*suppDistCorner*colDistCorner*cos(internalAngle)]));
            [SP_flip4, orders4, sing_flip] = symbolicStationaryPointsCorner(colDistCorner, suppDistCorner, fun, internalAngle, funSide, phaseCorner);
            real_sing_flip = mean(sing_flip);
            
            %construct singularity data
            R = @(yr) sqrt(colDistCorner^2 + (suppDistCorner + yr).^2 - 2*cos(internalAngle)*colPt.distSideL*(suppDistCorner + yr));
            logSingInfo_flip = singularity(real_sing_flip, Op.singularity, R);
            rectRad = .5*min(logSingInfo_flip.distFun(a_shift), logSingInfo_flip.distFun(b_shift));
            
            %now remove stationary points that are far away from
            %integration region:

            Xoscs = findNonOscBit(phaseCorner{1},0,width,kwave,minOscs);
            logSingInfo_flip_1.position = real(sing_flip(1));
            logSingInfo_flip_1.blowUpType='nearLog';
            logSingInfo_flip_1.distFun = @(r) abs(r-sing_flip(1));

            [x_, w_] = NonOsc45(0,Xoscs,kwave,Nquad,phaseCorner{1},logSingInfo_flip_1,Xoscs);
            I_1 = (w_.'*amp_corner(x_));
            if Xoscs>=width
                I_2 = 0;
            elseif ~isempty(SP_flip4)
                if SP_flip4<Xoscs
                    I_2 = PathFinderChebWrapGrad(logSingInfo_flip, Xoscs, width, kwave, Nquad, amp_corner, phaseCorner);
                else
                    if abs(sing_flip(1) - Xoscs) < dangerWidth
                        error('Singularity dangerously close and not being acknowledged');
                    end
                    if ~isempty(quadDataIn)
                        z = quadDataIn.z;
                        w = quadDataIn.w;
                    else
                        [ z, w ] = PathFinder( Xoscs, width, kwave, Nquad, phaseCorner,...
                                        'stationary points', SP_flip4, 'order', orders4, 'settlerad', ...
                                            rectRad, 'minOscs', minOscs, 'width', width);
                        quadDataOut.z = z;
                        quadDataOut.w = w;
                    end
                    I_2 = (w.'*amp_corner(z));
                end
            else
                I_2 = PathFinderChebWrapGrad(logSingInfo_flip,Xoscs,width,kwave,Nquad,amp_corner,phaseCorner);
            end
            I = I_1 + I_2;
    end
    
%    quadDataOut = struct('z',z,'w',w,'split',split,'z1a', z1a, 'w1a',w1a,'z1b', z1b, 'w1b', w1b);
            
end