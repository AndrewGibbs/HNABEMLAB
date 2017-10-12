function output = QuadTestCol( kwave,pMax, QuadIn )
%tests collocation matrices computed using different quadrature

display('Computing same discrete system with two different quadrature rules:');
fprintf('Brute force composite Gauss V %s\n',class(QuadIn));
    try

        %define the screen
        vertices=[0 0       %first vertex
                  1 0];     %second vertex

        %create 'edge' object for the screen
        Gamma=edge(vertices);

        %inident plane wave
        uinc=planeWave(kwave,[1 1]./sqrt(2));

        %with RHS data
        f=DirichletData(uinc,Gamma);

        %make an HNA basis on Gamma
        nLayers=2*(pMax+1); sigmaGrad=0.15; alphaDist=2;
        VHNA=HNAsingleMesh(Gamma,pMax,kwave,alphaDist, nLayers, sigmaGrad);
        DOFs=length(VHNA.el);
        %define the single layer 'operator' object
        S=SingleLayer(kwave,Gamma);

        %construct Geometrical optics approximation on Gamma
        GOA=GeometricalOpticsApprox(uinc,Gamma);

        X = getColPoints( VHNA );
        %X = ChebyshevRoots( DOFs, 'Tn', [Gamma.supp(1) Gamma.supp(2)] ).';

        preColMatrix=[];
        elCount=1;

        %create collocation matrix
        for b=VHNA.el
            Sb=S*b;
            preColMatrix=[preColMatrix Sb.col(X)];
            elCount=elCount+1;
        end

        %create RHSa, no integrals here
        %construct Geometrical optics approximation on Gamma
        GOA=GeometricalOpticsApprox(uinc,Gamma);% so can be computed exactly

        %create RHSb, these integrals will also have to be solved
        Sgoa=S*GOA;
        preColRHSb=Sgoa.col(X);

        % initialise brute force integration solver to test against
        CG=CompGaussBasic(kwave,5000,5000);

        ColMatrixCG=zeros(size(preColMatrix));
        ColMatrixNSD=zeros(size(preColMatrix));
        ColRHSb=zeros(DOFs,1); ColRHSbGC=ColRHSb; ColRHSbNSD=ColRHSb;
        for n=1:DOFs
            CG_=CG;
            QI_=QuadIn;
            fprintf('Computing column %d of %d\n',n,DOFs);
            for m=1:DOFs
                ColMatrixCG(m,n)=CG_.eval(preColMatrix(m,n));
                ColMatrixNSD(m,n)=QI_.eval(preColMatrix(m,n));
            end
            ColRHSbGC(n)=CG_.eval(preColRHSb(n));
            ColRHSbNSD(n)=QI_.eval(preColRHSb(n));
        end
        errL=MatMax(abs(ColMatrixCG-ColMatrixNSD));
        errR=MatMax(abs(ColRHSbGC-ColRHSbNSD));
        output=char('Matrix error: ', errL, 'RHS error: ', errR);
        display('Comparison complete');
    catch MatlabError
        display('An error occured');
        output=char('Matlab errored: ',MatlabError.identifier);
    end
end

