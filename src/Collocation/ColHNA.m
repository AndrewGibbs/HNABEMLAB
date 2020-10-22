function [v_N, GOA, colMatrix, colRHS, solveTime] = ColHNA(Operator, Vbasis, uinc, Gamma, varargin)
%computes oversampled collocation projection using HNA basis/frame

    %surpress warnings about cleared data in parfor loop:
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
    
    if iscell(uinc)
        squiggly = true;
    else
       squiggly = false;
       uinc_ = uinc;
       clear uinc;
       uinc{1} = uinc_;
    end
    num_incs = length(uinc);
    
    % -----------------------
    %defaults:
    Nquad = 15;
    colType='C';
    overSamplesPerMeshEl=1;
    messageFlag=false;
    standardBEMflag = false;
    standardQuadFlag = false;
    truncParam = 1e-8;
    weighting = false;
    symmetry_trick = false;
    % -----------------------
    
    for j=1:length(varargin)
        if ischar(varargin{j})
           lowerCaseArg=lower(varargin{j});
           switch lowerCaseArg
               case 'oversample'
                   overSamplesPerMeshEl=varargin{j+1};
               case 'coltype'
                   colType=varargin{j+1};
               case 'quadpoints'
                   Nquad=varargin{j+1};
               case 'max stationary point order'
                   maxSPorder=varargin{j+1};
               case 'progress'
                   messageFlag = true;
               case 'cg'
                   standardQuadFlag = true;
               case 'weight'
                   weighting = true;
               case 'trunc'
                   truncParam = varargin{j+1};
               case 'symmetry'
                   symmetry_trick = true;
                   
           end
        end
    end
    
    if isa(Vbasis,'hpStandardBasis')
       standardBEMflag = true;
    end
    
    %start timer
    tic;
    
    %get collocation points.
    Xstruct = getColPointsV2( Vbasis, overSamplesPerMeshEl, colType);
    %make a copy for RHSs, with wider domains
    Ystruct = Xstruct;
    for m = 1:length(Xstruct)
       Ystruct(m).distMeshL = Ystruct(m).distSideL;
       Ystruct(m).distMeshR = Ystruct(m).distSideR; 
    end
    
    %initialise main bits:
    colMatrix=zeros(length(Xstruct),length(Vbasis.el));
    numColPts = length(Xstruct);
    numBasEls = length(Vbasis.el);
    
    %symmetrySearch %not recommended for fractals with complex symmetry structure
    if symmetry_trick
        [colSymIndices,basisSymIndices] = getSymmetryIndices(Operator,numBasEls,numColPts);
    else
        colSymIndices = nan(numColPts,numBasEls);
    end
    
    if messageFlag
       fprintf('\nConstructing BEM matrix:');
    end
    parfor m=1:numColPts %can be parfor
        %fprintf('\nm');
        VbasisCopy = Vbasis;
        colMatrixCol = zeros(1,numBasEls);
        for n=1:numBasEls
            if isnan(colSymIndices(m,n))
            %manually do first entry of row
               if n==1
                   [colMatrixCol(n), quadData] = LHSquad(Operator, VbasisCopy.el(1), VbasisCopy.elEdge(1), Xstruct(m), Nquad,[]);
               elseif VbasisCopy.el(n).pm == VbasisCopy.el(n-1).pm && isequal(VbasisCopy.el(n).meshEl,VbasisCopy.el(n-1).meshEl)
                       % && VbasisCopy.el(n).gradIndex==VbasisCopy.el(n-1).gradIndex%added this line as supports can seem equal when they're not
                   %reuse quadrature from previous iteration of this loop,
                   %(phase and domain are the same)
                   colMatrixCol(n) = LHSquad(Operator, VbasisCopy.el(n), VbasisCopy.elEdge(n), Xstruct(m), Nquad, quadData);
               else
                   %get fresh quadrature data
                   [colMatrixCol(n), quadData] = LHSquad(Operator, VbasisCopy.el(n), VbasisCopy.elEdge(n), Xstruct(m), Nquad, []);
               end
               
               if Xstruct(m).side == VbasisCopy.elEdge(n)
                   colMatrixCol(n) = colMatrixCol(n) + Operator.Id(VbasisCopy.el(n), Xstruct(m).x);
               end
            end
        end
         colMatrix(m,:) = colMatrixCol;
        if messageFlag
            fprintf('\n\t%d/%d%',m,numColPts);
        end 
    end
    if messageFlag
       fprintf('\nConstructing RHS vector(s):');
    end
    %new seperate RHS bit, which allows multiple RHSs
    colRHS=zeros(length(Xstruct),num_incs);
    for ui_=1:num_incs
        colRHS_=zeros(length(Xstruct),1);
        f = Operator.get_RHS_data(uinc{ui_});
        %construct Geometrical optics approximation on Gamma
        GOA{ui_}=GeometricalOpticsApprox(uinc{ui_},Gamma);
        parfor m=1:numColPts
            fCopy = f;
            %with RHS data
            %colMatrix(m,:) = colMatrixCol;
            fX = fCopy.eval(Xstruct(m).x,Xstruct(m).side);
            if ~standardBEMflag
               OpPsiX = 0;
               for GOAedge = GOA{ui_}.suppEdges
                   OpPsiX = OpPsiX + RHSquad(Operator, GOA{ui_}.edgeComponent(GOAedge), GOAedge, Ystruct(m), Nquad);
               end
               colRHS_(m)  = fX - (OpPsiX + Operator.Id(GOA{ui_}.edgeComponent(Xstruct(m).side),Xstruct(m).x));
            else
               colRHS_(m)  = fX;
            end
            if messageFlag
                fprintf('\n\t%d/%d%',m,numColPts);
            end 
        end
        colRHS(:,ui_) = colRHS_;
    end
    
    %now fill in any gaps in matrix which can be done by exploiting
    %symmetry in structure
    for m=1:numColPts
       for n = 1:numBasEls
           if ~isnan(colSymIndices(m,n))
               colMatrix(m,n) = colMatrix(basisSymIndices(m,n),colSymIndices(m,n));
           end
       end
    end
    
    if weighting
        for j=1:length(Xstruct)
            w(j) = Xstruct(j).weight;
        end
        weightrix = (diag(sqrt(w)));
        %now weight LS system by weighted matrix
        colRHS = weightrix * colRHS;
        colMatrix = weightrix * colMatrix;
    end
    
    %use least squares with SVD to get coefficients
    if isnan(truncParam)
        coeffs = colMatrix\colRHS;
    elseif strcmp(truncParam,'inv')
        coeffs = inv(colMatrix)*colRHS;
    else
        %or, use Daan's homemade SVD:
        coeffs = pseudo_backslash(colMatrix, colRHS, truncParam);
    end
    solveTime = toc;
    if squiggly
        for m=1:num_incs
            v_N{m} = ProjectionFunction(coeffs(:,m),Vbasis);
        end
    else
        v_N=ProjectionFunction(coeffs,Vbasis);
        GOA = GOA{1};
    end
    
    if messageFlag
        toc
    end
    
end