function [v_N, GOA, colMatrix, colRHS] = ColHNA(Operator, Vbasis, uinc, Gamma, varargin)
%computes oversampled collocation projection using HNA basis/frame

    %with RHS data
    f=DirichletData(uinc,Gamma);
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);
    
    DOFs=length(Vbasis.el);
    
    % -----------------------
    %defaults:
    Nquad = 20;
    scaler=1;
    colType='C';
    overSamplesPerMeshEl=1;
    messageFlag=false;
    standardBEMflag = false;
    standardQuadFlag = false;
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
    [X, Xside, Xdist2a, Xdist2b] = getColPoints( Vbasis, overSamplesPerMeshEl, scaler, colType);
    
    %** should eventually find a way to partition the basis into mesh
    %elements with + or - phase, so that quadrature points can be reused.
    
    %initialise main bits:
    colRHS=zeros(length(X),1);
    colMatrix=zeros(length(X),length(Vbasis.el));
    for m=1:length(X)
        %manually do first entry of row
        [colMatrix(m,1), quadData] = colEval(Operator, Vbasis.el(1), Vbasis.elSide(1), X(m), Xside(m), Xdist2a(m), Xdist2b(m), Nquad,[], standardQuadFlag);
        for n=2:length(Vbasis.el)
           if Vbasis.el(n).pm == Vbasis.el(n-1).pm && isequal(Vbasis.el(n).supp,Vbasis.el(n-1).supp) && 2+2==5
               %reuse quadrature from previous iteration of this loop,
               %(phase and domain are the same)
               colMatrix(m,n) = colEval(Operator, Vbasis.el(n), Vbasis.elSide(n), X(m), Xside(m), Xdist2a(m), Xdist2b(m), Nquad, quadData);
               1+1;
           else
               %get fresh quadrature data
               [colMatrix(m,n), quadData] = colEval(Operator, Vbasis.el(n), Vbasis.elSide(n), X(m), Xside(m), Xdist2a(m), Xdist2b(m), Nquad,[], standardQuadFlag);
                1+1;
           end
        end 
        fX = f.eval(X(m),Xside(m));
        if ~standardBEMflag
           SPsiX =  colEval(Operator, GOA,GOA.illumSides, X(m), Xside(m), X(m), X(1+length(X)-m), Nquad,[], standardQuadFlag);
           colRHS(m)  = fX - SPsiX;
        else
           colRHS(m)  = fX;
        end
        if messageFlag
            fprintf('\n%.1f%%',100*m/length(X));
        end
    end
    
    %use least squares with Matlab's built in SVD to get coefficients
    coeffs=colMatrix\colRHS;
    v_N=Projection(coeffs,Vbasis);
    
    if messageFlag
        toc
    end
    
end