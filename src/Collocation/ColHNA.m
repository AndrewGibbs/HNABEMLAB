function [v_N, GOA, colMatrix, colRHS, colMatrix2, colRHS2] = ColHNA(Operator, Vbasis, uinc, Gamma, varargin)
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
    [X, Xside] = getColPoints( Vbasis, overSamplesPerMeshEl, scaler, colType);
    
    %** should eventually find a way to partition the basis into mesh
    %elements with + or - phase, so that quadrature points can be reused.
    
    %initialise main bits:
    colRHS=zeros(length(X),1);
    colMatrix=zeros(length(X),length(Vbasis.el));
%     
%     %start the main double loop:
%     for m=1:length(X)
%         for n=1:length(Vbasis.el)
%            colMatrix(m,n) = colEval(Operator, Vbasis.el(n), Vbasis.elSide(n), X(m), Xside(m), Nquad);
%         end 
%         if ~standardBEMflag
%             colRHS(m) = f.eval(X(m),Xside(m)) - colEval(Operator,GOA,GOA.illumSides,X(m),Xside(m),Nquad);
%         else
%             colRHS(m) = f.eval(X(m),Xside(m));
%         end
%         if messageFlag
%             fprintf('\n%.1f%%',100*m/length(X));
%         end
%     end
    colMatrix2 = colMatrix;
    colRHS2 = colRHS;
    for m=1:length(X)
        for n=1:length(Vbasis.el)
           [colMatrix(m,n), colMatrix2(m,n) ] = colEval(Operator, Vbasis.el(n), Vbasis.elSide(n), X(m), Xside(m), Nquad);
        end 
        fX = f.eval(X(m),Xside(m));
        if ~standardBEMflag
           [    SPsiX, SPsiX2] =  colEval(Operator,GOA,GOA.illumSides,X(m),Xside(m),Nquad);
           colRHS(m)  = fX - SPsiX;
           %colRHS2(m) = fX - SPsiX2;
           colRHS2(m) = SPsiX2;
        else
           colRHS(m)  = fX;
           colRHS2(m) = fX;
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