function [v_N, VHNA, X, coeffs] = HNAColOversample( kwave,Gamma,uinc,pMax,extraPoints,scaler, messageFlag )
%function which returns oversampled collocation HNA BEM
    if nargin<=5
        scaler=1;
    end
    if nargin<=6
        messageFlag=0;
    end
    %with RHS data
    f=DirichletData(uinc,Gamma);

    %make an HNA basis on Gamma
    nLayers=2*(pMax+1); sigmaGrad=0.15;
    VHNA=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);
    DOFs=length(VHNA.el);
    %define the single layer 'operator' object
    S=SingleLayer(kwave,Gamma);

    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);

    X = getColPoints( VHNA, extraPoints, scaler);
    COLs=length(X);

    preColMatrix=[];
    elCount=1;

    %create collocation matrix
    for b=VHNA.el
        if messageFlag==1
           fprintf('Creating column %d of %d\n',elCount,DOFs);
        end
        Sb=S*b;
        preColMatrix=[preColMatrix Sb.col(X)];
        elCount=elCount+1;
    end

    %create RHSa, no integrals here
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);% so can be computed exactly
    RHSa=f.eval(X);

    %create RHSb, these integrals will also have to be solved
    Sgoa=S*GOA;
    preColRHSb=Sgoa.col(X);

    % initialise integration solver
    NSD=NSDlinearPhase(20,20);

    ColMatrixNSD=zeros(size(preColMatrix));
    ColRHSb=zeros(COLs,1);
    for m=1:COLs
        if messageFlag==1
           fprintf('Computing row %d of %d\n',m,COLs);
        end
        for n=1:DOFs
            ColMatrixNSD(m,n)=NSD.eval(preColMatrix(m,n));
        end
        ColRHSb(m)=NSD.eval(preColRHSb(m));
    end
   % fprintf('Condition number of matrix: %e',cond(ColMatrixNSD));
    ColRHS=RHSa-ColRHSb;
    warning('off','MATLAB:rankDeficientMatrix');
    coeffs=ColMatrixNSD\ColRHS;
    v_N=Projection(coeffs,VHNA);

end

