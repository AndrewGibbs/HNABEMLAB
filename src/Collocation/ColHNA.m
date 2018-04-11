function [v_N, GOA, colMatrix, colRHS] = ColHNA(Operator, HNAbasis, uinc, Gamma, varargin)
%computes oversampled collocation projection using HNA basis/frame
    kwave = uinc.kwave;
    %with RHS data
    f=DirichletData(uinc,Gamma);
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);
    
    DOFs=length(HNAbasis.el);
    
    %defaults:
    Nquad = 15;
    scaler=1;
    colType='C';
    overSamplesPerMeshEl=1;
    maxSPorder=0;
    
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
           end
        end
    end
    
    
    %get collocation points.
    X = getColPoints( HNAbasis, overSamplesPerMeshEl, scaler, colType);
    
    %** should eventually find a way to partition the basis into mesh
    %elements with + or - phase, so that quadrature points can be reused.
    
    %initialise main bits:
    colRHS=zeros(length(X),1);
    colMatrix=zeros(length(X),length(HNAbasis.el));
    
    %start the main double loop:
    for m=1:length(X)
        for n=1:length(HNAbasis.el)
           colMatrix(m,n) = colEval(Operator,HNAbasis.el(n),X(m),maxSPorder,Nquad);
        end 
        colRHS(m) = f.eval(X(m)) - colEval(Operator,GOA,X(m), maxSPorder,Nquad);
        fprintf('\n%.1f%%',100*m/length(X));
    end
    
    %use least squares with Matlab's built in SVD to get coefficients
    coeffs=colMatrix\colRHS;
    v_N=Projection(coeffs,HNAbasis);
    
end