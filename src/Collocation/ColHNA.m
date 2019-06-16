function [v_N, GOA, colMatrix, colRHS] = ColHNA(Operator, Vbasis, uinc, Gamma, varargin)
%computes oversampled collocation projection using HNA basis/frame

    %with RHS data
    f=DirichletData(uinc,Gamma);
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);
    
    DOFs=length(Vbasis.el);
    
    weighting = false;
    
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
    
    %start timer, if user is interested
    if messageFlag
        tic;
    end
    
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
    parfor m=1:numColPts %can be parfor
        %fprintf('\nm');
        VbasisCopy = Vbasis;
        fCopy = f;
        colMatrixCol = zeros(1,numBasEls);
        for n=1:numBasEls
        %manually do first entry of row
           if n==1
               [colMatrixCol(n), quadData] = colEvalV3(Operator, VbasisCopy.el(1), VbasisCopy.elEdge(1), Xstruct(m), Nquad,[], standardQuadFlag);
           elseif VbasisCopy.el(n).pm == VbasisCopy.el(n-1).pm && isequal(VbasisCopy.el(n).supp,VbasisCopy.el(n-1).supp)
               %reuse quadrature from previous iteration of this loop,
               %(phase and domain are the same)
               colMatrixCol(n) = colEvalV3(Operator, VbasisCopy.el(n), VbasisCopy.elEdge(n), Xstruct(m), Nquad, quadData, standardQuadFlag);
           else
               %get fresh quadrature data
               [colMatrixCol(n), quadData] = colEvalV3(Operator, VbasisCopy.el(n), VbasisCopy.elEdge(n), Xstruct(m), Nquad,[], standardQuadFlag);
           end
        end
        %integral(@(t) Operator.kernel(abs(t-Xstruct(m).x)).*VbasisCopy.el(n).eval(t), VbasisCopy.el(n).a, VbasisCopy.el(n).b);
        %abs(colMatrixCol(n)-colEvalV2(Operator, VbasisCopy.el(n),VbasisCopy.elEdge(n), Xstruct(m), Nquad,[],standardQuadFlag))>1e-8
        colMatrix(m,:) = colMatrixCol;
        fX = fCopy.eval(Xstruct(m).x,Xstruct(m).side);
        if ~standardBEMflag
           SPsiX = 0;
           for GOAedge = GOA.suppEdges
               SPsiX = SPsiX + colEvalV3(Operator, GOA.edgeComponent(GOAedge), GOAedge, Ystruct(m), Nquad,[], standardQuadFlag);
           end
           colRHS(m)  = fX - SPsiX;
        else
           colRHS(m)  = fX;
        end
        if messageFlag
            fprintf('\n%d/%d%',m,numColPts);
        end
    end
    
%     for m=1:numColPts
%         for n=1:numBasEls
%             colMatrixBF(m,n) = integral(@(t) Operator.kernel(abs(t-Xstruct(m).x)).*Vbasis.el(n).eval(t), Vbasis.el(n).a, Vbasis.el(n).b);
%         end
%     end
    
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
    v_N=ProjectionFunction(coeffs,Vbasis);
    
    if messageFlag
        toc
    end
    
end