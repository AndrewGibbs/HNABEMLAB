function [v_N, GOA, colMatrix, colRHS, solveTime] = ColHNA(Operator, Vbasis, uinc, Gamma, varargin)
%computes oversampled collocation projection using HNA basis/frame

    %surpress warnings about cleared data in parfor loop:
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
    
    %with RHS data
    f=DirichletData(uinc,Gamma);
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);
    
    DOFs=length(Vbasis.el);
    
    % -----------------------
    %defaults:
    Nquad = 15;
    scaler=1;
    colType='C';
    overSamplesPerMeshEl=1;
    messageFlag=false;
    standardBEMflag = false;
    standardQuadFlag = false;
    truncParam = 1e-8;
    symmetrySearch = false;
    weighting = false;
    symmetry_trick = true;
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
    colRHS=zeros(length(Xstruct),1);
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
    for m=1:numColPts %can be parfor
        %fprintf('\nm');
        VbasisCopy = Vbasis;
        fCopy = f;
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
            end
        end
        colMatrix(m,:) = colMatrixCol;
        fX = fCopy.eval(Xstruct(m).x,Xstruct(m).side);
        if ~standardBEMflag
           SPsiX = 0;
           for GOAedge = GOA.suppEdges
               SPsiX = SPsiX + RHSquad(Operator, GOA.edgeComponent(GOAedge), GOAedge, Ystruct(m), Nquad);
           end
           colRHS(m)  = fX - SPsiX;
        else
           colRHS(m)  = fX;
        end
        if messageFlag
            fprintf('\n\t%d/%d%',m,numColPts);
        end
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
    
    %use least squares with Matlab's built in SVD to get coefficients
    if isnan(truncParam)
        coeffs = colMatrix\colRHS;
    elseif strcmp(truncParam,'inv')
        coeffs = inv(colMatrix)*colRHS;
    else
        %or, use Daan's homemade SVD:
        coeffs = pseudo_backslash(colMatrix, colRHS, truncParam);
    end
    solveTime = toc;
    v_N=ProjectionFunction(coeffs,Vbasis);
    
    if messageFlag
        toc
    end
    
end