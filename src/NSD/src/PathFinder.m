function [ X, W, pathData ] = PathFinder( a,b,freq,N,G,varargin)
%returns O(N) weights w and nodes z which evaluate an oscillatory integral over the interval [a,b], with
%frequency parameter k, and phase g(z). The input G is an array of
%anonymous functions, such that G{n} corresponds to g^(n-1)(z). For
%non-degenerate stationary points, we only require as high as n=3;
%stationary points of higher order require more derivatives.

%Where important information about the phase g(x) is unknown, this code
%will approximate the information using and ODE45.

%if there is a stationary point just outside of [a,b] on the real line what do we do about
%it? It will get picked up by the root finder. If we can plot the level
%curvers from h_a(Pmax) to h_b(Pmax), can check for stationary points
%inside of this region... maybz

%if the points on the SD path end up outside of the initial rectangle, then
%there might be stationary points nearby which should be considered; in
%fact it may be sensible to start over with a larger rectangle, to be
%totally robust.

%add all of the necessary code:
%StandardPaths();

%known bugs:
% if a singularity has the same position as a stationary point, PathFinder
% doesn't return a complete path

    %% ------------------------------------------------------------------%
    % -------------------------- KEY PARAMETERS ------------------------ %
    %--------------------------------------------------------------------%
    
    %loads of different values for 'almost zero'
        %tolerance for ODE45 solver:
        RelTol=1E-12;
        %tolerance for two statoionary points to become clumped together:
        RectTol=min(10000*(1/freq),1e-12);
        TaylorRad = 1E-4;%2*max(1E-12,RectTol);
        %distance from a singularity to be considered nearly singular:
        nearSingThresh = 0.2;%0.1;
        %distance from an integer to be considered not an integer
        intThresh=0.01;
        %number of oscillations to actually be considered oscillatory:
        minOscs = 2;
        tweakThresh = eps*10;
        tweakableJacobianThresh = 2;
        %tweakPaths = false;
        
    %number of quadrature points used for the integrals inside of this code
        RectQuadPts = 50;
    %largest number of quad points to use
        quadMax = 50;
        
    %tests for nearly finite paths on or off:
        nearlyFiniteTests=false;
        
    %flag for plotting stuff as we go
        visuals=false;
        gAnalytic=true;
        
    %default rectangle radius - dependence on frequency not totally clear
    %yet
        rectRad=.5/freq;
        
    %default f stuff:
        fSingularities=[]; fSPorders=[];
        gSingularities = [];
        
    %default settle radius
        settleRad=[];
        
        %default inf flags
        ainf=false;
        binf=false;
        
        width = b-a;
        
        getPathData = false;
        widthIndex = [];
    
    %% ------------------------------------------------------------------%
    % -------------------- INTERPRET USER INPUT ------------------------ %
    %--------------------------------------------------------------------%

    %check that G is full of function handles - a misplaced space can mess
    %this up
    for j=1:length(G)
        if ~isa(G{j}, 'function_handle') 
            error(sprintf('%dth entry of 5th input argument not function handle',j));
        end
    end
    
    if length(varargin)==1 && ~ischar(varargin{1})
        %glitchy Matlab varagin thing, only an issue when this function 
        %call's it's self recursively. Fix it here:
        varargin=varargin{1};
        %MESSES UP WHEN YOU SEND MATLAB ONLY ONE OPTION, DUH
    end
    %use NaN to store things which are not yet defined, as opposed to
    %empty (which is denoted [], and may be specified by the user):
    gStationaryPoints=NaN;
    fSingularities=[]; fSingularitiesObj=[]; %singularities in f, not g
    Mf=1; fTest=false;
    interiorSingularities=[]; fSingularities=[];
    %check through optional inputs
    for j=1:length(varargin)
        if ischar(varargin{j})
           lowerCaseArg=lower(varargin{j});
           switch  lowerCaseArg
               case 'stationary points'
                   gStationaryPoints=varargin{j+1};
                   if isempty(gStationaryPoints)
                       gSPorders=[];
                   end
                   %RectTol = 1E-16; %this can be much smaller now, given that stationary point values are exact
               case 'pathData'
                   getPathData = true;
               case 'fsingularities'
                   fSingularitiesObj=varargin{j+1};
                    %check for singularities in (real) interior of integration range:
                    for s=fSingularitiesObj
                        if imag(s.position)==0
                           if  a<s.position && s.position<b
                               interiorSingularities=[interiorSingularities s.position];
                           end
                        else %keep hold of these for later, incase deformation crosses them
                            fSingularities=[fSingularities s.position];
                        end
                    end
               case 'gSingularities'
                   gSingularities=varargin{j+1};
               case 'order'
                   gSPorders=varargin{j+1};               
%                case 'mf' %upper bound of non-oscillatory f in complex plane
%                    Mf=varargin{j+1};
% %                case 'ginv' %user has provided inverse of g(x)
% % 
% %                case 'gderinv' %user has provided inverse of g'(x)
% %                    
               case 'ganalytic'
                   gAnalytic=varargin{j+1};
               case 'visuals on'
                   visuals=true;
                   hold on;
               case 'ftest' %not necessary, but user can input f, to find singularities etc
                   fTest=true;
                   F=varargin{j+1};
                   if length(F) <2
                       error('need f(x) and its derivative, in cell form');
                   end
               case 'rectrad'
                   rectRad=varargin{j+1};
               case 'settlerad' %radius outside of which function settles down
                   settleRad = varargin{j+1};
                   rectRad = settleRad;
               case 'ainf'
                   ainf=true;
                   width=inf;
               case 'binf'
                   binf=true;
                   width=inf;
               case 'gpolycoeffs'
                   polyCoeffs = varargin{j+1};
                   [G, gStationaryPoints, gSPorders] = NSDeetsFromPoly(polyCoeffs, RectTol);
               case 'width'
                   width = varargin{j+1};
                   widthIndex = j+1;
               case 'minoscs'
                   minOscs = varargin{j+1};
               case 'TaylorRad'
                   TaylorRad = varargin{j+1};
               case 'RectTol'
                   RectTol = varargin{j+1};
               case 'nearSingThresh'
                   RectTol = varargin{j+1};
           end
       end
    end
    if width <=0
        error('Integration interval must have positive width');
    end
    if isinf(minOscs)
        if isinf(width)
            error('Cannot brute force an infinite contour integral, minOscs must be finite if ainf or binf options are turned on');
        end
        [X, W] = NonOsc45(a,b,freq,N,G{1},fSingularitiesObj, width);
        return;
    end
    
    %% ------------------------------------------------------------------%
    % -------------------------- SINGULARITIES ------------------------- %
    %--------------------------------------------------------------------%
         
    %scan for singularities, if requested:
    if fTest
        [~, fSingularities, fSPorders, ~] = findZerosSingsRect( F{1}, F{2}, initRect, RectTol, RectQuadPts , visuals);
    end
    
    %if there are any interior singularities, run iteratively with them at
    %the endpoint in each case
    if ~isequal([],interiorSingularities) 
        %singularities in interior, so let's make like a banana and split
        intervalSplit=[a interiorSingularities b];
        X=[]; W=[];
        for j=1:(length(intervalSplit)-1)
            %call self recursively, without interior singularities:
            %vararginCopy=varargin;
            if isempty(widthIndex)
                [ X_, W_ ] = PathFinder( intervalSplit(j),intervalSplit(j+1),freq,N,G,varargin);
            elseif length(width) ~= (length(intervalSplit)-1)
                %warning('User has provided width of integration interval, however this cannot be used, as interval has been split at singularity');
                varargin_ = varargin;
                varargin_{widthIndex-1} = 'nothing useful'; %this will avoid the width value getting used, and it'll just get ignored instead
                [ X_, W_ ] = PathFinder( intervalSplit(j),intervalSplit(j+1),freq,N,G,varargin_);
            else
                varargin_ = varargin;
                varargin_{widthIndex} = width(j);
                [ X_, W_ ] = PathFinder( intervalSplit(j),intervalSplit(j+1),freq,N,G,varargin_);
            end
            X=[X; X_;]; W=[W; W_;];
        end
        %and end routine 'early'
        return;
    end 
    
    
    %% ------------------------------------------------------------------- 
    % ---------------- STAIONARY and CCRITICAL POINTS ------------------ %
    %--------------------------------------------------------------------%
    
    if isnan(gStationaryPoints) %no stationary points specified by user
        [gStationaryPoints, gSingularities, gSPorders,~] = getStationaryPoints(a,b,rectRad,...
                                                            gAnalytic, G, RectTol, RectQuadPts , visuals,...
                                                            settleRad, intThresh);
    else %user has provided stationary points, but now let's throw away those which are nowhere near [a,b]
        %[gStationaryPoints, gSPorders] = stationaryPointsCheck(G,gStationaryPoints,gSPorders);
        [gStationaryPoints, gSPorders] = pruneStationaryPoints(a,b,rectRad,gStationaryPoints, gSPorders, ainf, binf, settleRad);
    end
    
    %locate branch points ** this probably needs either deleting or
    %updating?
    branchPoints = getIntegrandBranchPoints(gSPorders, fSPorders, intThresh);
    
    %combine stationary points with [a,b] to make 'critical points'
    [criticalPoints, pathPowers] = makeCriticalPoints(a,b,gStationaryPoints,gSPorders+1, RectTol, gAnalytic, ainf, binf);
    
    
    %% Check that function is actually oscillatory:
    oscs = max(oscCount( a,b,freq,gStationaryPoints,G{1} ));
    if max(oscs) < minOscs && ~isinf(width)
        [X, W] = NonOsc45(a,b,freq,N,G{1},fSingularitiesObj, width,oscs);
        return;
    end
    
    %% ------------------------------------------------------------------%
    % -------------------------- COMPUTE PATHS ------------------------- %
    %--------------------------------------------------------------------%
    
    % if not enough derivatives provided:
    if length(G)< (round(max(pathPowers))+1)
        error('Need at least %d derivatives of g(x) to compute paths',round(max(pathPowers))+1);
    end
    
    %------now do a load of bodging for the finite path stuff------%
%     
    % **** currently assume that there are no singularities on finite paths
    
    % ***** current code does not account for finite paths stemming from
    % different branches of same stationary point. Will need to replace the
    % vector of distances by a matrix, add an extra loop in
    % getBranchesOfFinitePaths, and also account for the cases where one
    % path goes through three stationary points... aaghsifdhdheh@##$%
    
    % ***** not sure if current finite path code can handle finite paths
    % with multiplicity one at one (or both) ends
    
     [prePathLengthsVec, pathEndIndex, ~, ~] = finitePathTest( criticalPoints, G , pathPowers);
     
     [hFinite, dhdpFinite, Wfinite, FPindices] = getBranchesOfFinitePaths(prePathLengthsVec, pathEndIndex, criticalPoints, pathPowers, G, freq, N, [], RelTol);
    
    %--- end of finite path bodging ---------------------------------%
    
    %now loop over all paths:
    fullIndex=0;    %this will be a unique combination of critPointIndex and branchIndex
    for critPointIndex=1:length(criticalPoints)
        
        %set ICs for SD path differential equation, which will determine SD path:
        ICs = NSDpathICv3( pathPowers(critPointIndex), G, criticalPoints(critPointIndex));
        %initialise this vector map, n'th entry is all the paths indices
        %which start at the n'th critical point
        CritPointToPathInds{critPointIndex}=[];
        
        for branchIndex=1:pathPowers(critPointIndex)
            fullIndex=fullIndex+1;
        
            if ~FPindices{critPointIndex}(branchIndex)
% 
%                 COVsingularities=fSingularitiesObj;
%                 for s=1:length(fSingularitiesObj)
%                     COVsingularities(s).position=(fSingularitiesObj(s).position-criticalPoints(critPointIndex))*(freq^(1/pathPowers(critPointIndex)))/1i;
%                     %THIS MAP only works for real singularities
%                 end

                %get weights and nodes
                    hNearSing = nearlySingularPhaseInverse(criticalPoints(critPointIndex), gStationaryPoints, nearSingThresh); 
                    [x{critPointIndex, branchIndex}, w{critPointIndex, branchIndex}]=pathQuadV2( criticalPoints(critPointIndex), pathPowers(critPointIndex), [fSingularitiesObj hNearSing], min(quadMax,N*pathPowers(critPointIndex)), freq , nearSingThresh);
                    P0=[0; (x{critPointIndex, branchIndex}./(freq^(1/pathPowers(critPointIndex))))];
                    %**** in above case, weights and nodes are independent of
                    %branchIndex
                
                %make Taylor series approximation close to stationary
                %point:
                hTaylor{critPointIndex,branchIndex} = TaylorPath(pathPowers(critPointIndex), G, criticalPoints(critPointIndex), branchIndex, false);
                %now solve IVP:
                [~,H] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(critPointIndex)-1,G, ICs{branchIndex}, false, TaylorRad, hTaylor{critPointIndex,branchIndex}), P0, ICs{branchIndex}, odeset('RelTol', RelTol) );
                
                
                %NEED TO INPUT HIGHER ORDER DE's INTO NSDpathODE.m
                %H(:,1) contains h'(p), H(:,2) contains h(p), will throw away initial
                %points

                if pathPowers(critPointIndex)>1
                    % h'(p) given as output of ODE45, so use it
                    h{critPointIndex,branchIndex}=H(2:end,1);
                    dhdp{critPointIndex,branchIndex}=H(2:end,2);
                else
                    %re-insert h(p) into DE for h'(p)*
                    h{critPointIndex,branchIndex}=H(2:end);
                    dhdp{critPointIndex,branchIndex}=1i./G{2}(H(2:end)); %back into the ODE 
                end
                
                if max(abs(dhdp{critPointIndex,branchIndex}))> tweakableJacobianThresh
                    tweakPaths = true;
                else
                    tweakPaths = false;
                end
                
                if tweakPaths
                    [h{critPointIndex,branchIndex}, dhdp{critPointIndex,branchIndex}]= pathAdjust(G, P0(2:end), h{critPointIndex,branchIndex}, criticalPoints(critPointIndex), pathPowers(critPointIndex), tweakThresh);
                end
                
                %display(max(abs(dhdp{critPointIndex,branchIndex})));

                %check that ODE45 gave a path that actually decays
                %exponentially:
%                 [divTF,nFinal] = divergenceTest(h{critPointIndex,branchIndex},G{1});
%                 if false%divTF
%                     h{critPointIndex,branchIndex}=h{critPointIndex,branchIndex}(1:nFinal);
%                     dhdp{critPointIndex,branchIndex}(1:nFinal);
%                 end
                
                %now check if the path we've just created is 'nearly
                %finite'
                if nearlyFiniteTests
                    [nearlyFinite, nfParamPoint, nfCritPointIndex] = nearlyFiniteCheck(P0, h{critPointIndex,branchIndex}, dhdp{critPointIndex,branchIndex}, criticalPoints(critPointIndex), criticalPoints, freq);
                else
                    nearlyFinite=false;
                end
                
                if nearlyFinite
                    %mash this nearly finite path into a finite path.
                    % update the FPindices vector
                    FPindices{critPointIndex}(branchIndex)=nfCritPointIndex;
                    [hFinite{critPointIndex,branchIndex}, dhdpFinite{critPointIndex,branchIndex}, Wfinite{critPointIndex,branchIndex}]...
                            = nearlyFinitePathFix(nfParamPoint, criticalPoints(nfCritPointIndex), criticalPoints(critPointIndex), pathPowers(critPointIndex),...
                                ICs{branchIndex}, G, freq, N, COVsingularities, RelTol);
                    X_{fullIndex} = hFinite{critPointIndex,branchIndex};
                    W_{fullIndex} = dhdpFinite{critPointIndex,branchIndex}.*Wfinite{critPointIndex,branchIndex};
                    %
                else
                    %this adjusts for small deviations from the path:
                    weightWatchers = exp(1i*freq*(G{1}(h{critPointIndex,branchIndex})-1i*P0(2:end).^pathPowers(critPointIndex) - G{1}(criticalPoints(critPointIndex))));
                    %absorb h'(p) and other constants into weights.
                    W_{fullIndex}=...
                        (1/(freq^(1/pathPowers(critPointIndex))))*exp(1i*freq*G{1}(criticalPoints(critPointIndex)))...
                        .*dhdp{critPointIndex,branchIndex}.*w{critPointIndex,branchIndex}...
                        .*weightWatchers;

                    X_{fullIndex}=h{critPointIndex,branchIndex}; 
                end
            else %finite paths have already been computed, necessarily to detect if they were finite
                %finite paths shouldn't diverge (in IVP sense) - else they become infinite
                %paths anyway
                X_{fullIndex}=hFinite{critPointIndex,branchIndex};
                W_{fullIndex}=dhdpFinite{critPointIndex,branchIndex}.*Wfinite{critPointIndex,branchIndex}*exp(1i*freq*G{1}(criticalPoints(critPointIndex)));
                %
            end
            if sum(isnan(X_{fullIndex}))>0
                error('NaNs in the complex plane. Possible explanation - a hidden stationary point?');
            end
            P(fullIndex,2)=X_{fullIndex}(end);   
            P(fullIndex,1)=criticalPoints(critPointIndex);
            CritPointToPathInds{critPointIndex}=[CritPointToPathInds{critPointIndex} fullIndex];
            PathIndsToCritPoint{fullIndex}=[critPointIndex branchIndex];
            CritPointEndPoint{critPointIndex}(branchIndex)=X_{fullIndex}(end);
        end
                
    end
    
    numPathsSD=fullIndex; %record total number of paths
    
    %% there is a difference between knowing the path, and walking the path%
    
    %now make a useful array which is indexed by the fullIndex, and returns
    %the fullIndex of every connected finite path
    fullIndex=0;
    for critPointIndex=1:length(criticalPoints)
        for branchIndex=1:pathPowers(critPointIndex)
            fullIndex=fullIndex+1;
            if FPindices{critPointIndex}(branchIndex)>0
                FPfullIndices{fullIndex}=CritPointToPathInds{FPindices{critPointIndex}(branchIndex)};
            else
                FPfullIndices{fullIndex}=[];
            end
        end
    end
    
    [X, W, finalPath,finalPathSign] = choosePathV2(a,b,criticalPoints, CritPointEndPoint, FPindices, P, G, freq, RectQuadPts, numPathsSD, visuals, X_, W_, hFinite, settleRad, ainf, binf,1);
    
    if getPathData
        %output more information so that data can be used as a starting
        %point for iterative procedure determining pertubations in phase
        for n = finalPath
            %get old coordinates back
            critPointIndex=PathIndsToCritPoint{n}(1);
            branchIndex=PathIndsToCritPoint{n}(2);
            %now get data which may be recycled
            pathData{n}.complexNodes = X_{n};
            pathData{n}.GaussNodes = x{critPointIndex, branchIndex};
            pathData{n}.GaussWeights = w{critPointIndex, branchIndex};
            pathData{n}.criticalPoint = criticalPoints(n);
            pathData{n}.order = gSPorders(n);
            pathData{n}.inOut = finalPathSign(n);
            if FPindices{critPointIndex}(branchIndex)
                pathData{n}.finitePath = true;
            else
                pathData{n}.finitePath = false;
            end
        end
    else
        pathData = [];
    end
    
    if ~gAnalytic || ~isempty(fSingularities)
        %add tiny circles around singular points if they lie inside the region
        %of deformation:
        [XsmallDisk, WsmallDisk] = handleSingularities(fSingularities, gSingularities, branchPoints, freq, G, N, visuals);
        X=[X; XsmallDisk];  W=[W; WsmallDisk];
    end
    
    %stop adding stuff to the lovely diagram
    if visuals
        hold off;
    end
end