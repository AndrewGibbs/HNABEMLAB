function [I, quadDataOut] = colEvalLinear(Op,fun, funSide, colPt, Nquad, quadDataIn)
%linearSpecialCase = false;
%function which evalutes integral Sf(x), essentially just a wrapper for
%NSD45 which ensures that the phase is always the analytic continuation of
%|x(s)-y(t)|
    
    if nargin <= 6
        CGflag = false;
    end
        
        kwave = Op.kwave;
        
        %get endpoints of support of function
        supp = fun.getSupp(funSide);
        a = supp(1);
        b = supp(2);
      
        %analytic extension of non-osc components of kernel:
        amp_a = @(y) Op.kernelNonOscAnal(colPt.x, y, true, colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_b = @(y) Op.kernelNonOscAnal(colPt.x, y, false, colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide);
        amp_a_flip = @(r) Op.kernelNonOscAnal(r, 0, true, colPt.side, funSide).* fun.evalNonOscAnal(colPt.x - r, funSide);
        amp_b_flip = @(r) Op.kernelNonOscAnal(r, 0, true, colPt.side, funSide).* fun.evalNonOscAnal(colPt.x + r, funSide);
        
        %now the more general amp, for when there is no branch in [a,b]
        
        maxSPorder = 0;
        amp = @(y) Op.kernelNonOscAnal(colPt.x, y, [], colPt.side, funSide) .* fun.evalNonOscAnal(y, funSide);
        phase = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, [], maxSPorder+1);
        phase_a = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, true, maxSPorder+1);
        phase_b = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, false, maxSPorder+1);
        %phase_b_flip = OpFunAddPhase(Op, fun, funSide, colPt.x, colPt.side, false, maxSPorder+1);
        
        for n = 1:length(phase_b)
            phase_a_flip{n} = @(r) (-1)^(n+1)*phase_a{n}(colPt.x-r);
            phase_b_flip{n} = @(r) phase_b{n}(r+colPt.x);
        end
        
        
        if length(fun.domain.L )>1
            L = fun.domain.L(funSide);
        else
            L = fun.domain.L;
        end
            
    if funSide == colPt.side
        %same side singularity:
        distFun = @(t) abs(colPt.x - t);
        logSingInfo=singularity(colPt.x, Op.singularity, distFun);
        
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
            end

           if ~isempty(quadDataIn)
                w1a = quadDataIn.w1a;
                w1b = quadDataIn.w1b;
                z1a = quadDataIn.z1a;
                z1b = quadDataIn.z1b;
                quadDataOut = quadDataIn;
            else
                logSingInfo_flip_a = logSingInfo;
                logSingInfo_flip_a.position = 0;
                logSingInfo_flip_a.distFun = @(r) abs(r);
                alphaPhase_a = phase_a_flip{2}(1);
                missingPhaseShift_a = phase_a_flip{1}(1) - alphaPhase_a;
                [ z1a, w1a ] = LinearNSD(colPt.distMeshL, kwave, alphaPhase_a, 0, Nquad);
                w1a = w1a*exp(1i*kwave*missingPhaseShift_a);

                logSingInfo_flip_b = logSingInfo;
                logSingInfo_flip_b.position = 0;
                logSingInfo_flip_b.distFun = @(r) abs(r);
                
                alphaPhase_b = phase_b_flip{2}(1);
                missingPhaseShift_b = phase_b_flip{1}(1) - alphaPhase_b;
                [ z1b, w1b ] = LinearNSD(colPt.distMeshR, kwave, alphaPhase_b, 0, Nquad);
                w1b = w1b*exp(1i*kwave*missingPhaseShift_b);
                
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
                            xi = (fun.meshEl.distR + fun.meshEl.width) - colPt.distSideR;
                        end
                    end
                    
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
                            missingPhaseShift = phase_star{1}(0) + a_star*alphaPhase;
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
            end
                    
            quadDataOut.w_ = w_;
            quadDataOut.z_ = z_;
                
            %and evaluate integral:
            I = w_.'*amp_star(z_);
            %now store in correct form:
%             if colPt.x <= a
%                 w1b = w_;
%                 z1b = z_;
%             else
%                 w1a = w_;
%                 z1a = z_;
%             end
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
                
            %need to fix rounding errors in (e.g.) "betaConst-alphaPhase*b",
                    % when these are all large, get different values to
                    % PathFinder and integral(.,.,.). However, the abs
                    % value of the weights is the same, suggesting that
                    % the error is only a rotation in the complex
                    % plane.
                alphaPhase = phase_shift{2}(0);
                betaConst = phase_shift{1}(0);s

                if alphaPhase>0
                    [ z_, w_ ] = LinearNSD(fun.suppWidth, kwave, alphaPhase, shiftedCol, Nquad);
                    w_ = w_*exp(1i*kwave*betaConst)*exp(1i*kwave*a*alphaPhase);
                else
                    [ z_, w_ ] = LinearNSD(fun.suppWidth, kwave, -alphaPhase, -shiftedCol, Nquad);
                    w_ = w_*exp(1i*kwave*(betaConst-alphaPhase*b));
                    z_ = b-z_;
                end
            end
        end
        quadDataOut.w_ = w_;
        quadDataOut.z_ = z_;
        
        %and evaluate integral:
        I = w_.'*amp_shift(z_);
            
    else %no branch in phase
        
         I = integral(@(z) amp(z).*exp(1i*kwave*phase{1}(z)), a, b, 'arrayValued', true,'RelTol',1e-13);
         quadDataOut = [];
            
end