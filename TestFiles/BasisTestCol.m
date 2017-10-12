function  resultsData  = BasisTestCol( kwave,pMax )

display('Solving same problem with different bases, to compare output');
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
    
    nLayers=2*(pMax+1); sigmaGrad=0.15; alphaDist=2;
    %make an HNA basis on a single mesh on Gamma
    VHNAs=HNAsingleMesh(Gamma,pMax,kwave,alphaDist, nLayers, sigmaGrad);
    %now make an overlapping mesh HNA basis
    VHNAo=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);
    %and a high order standard-hp-BEM basis to compare against
    hMax=2*pi/(2*kwave);
    Vhp=hpStandardBasis(Gamma, 8, hMax, 16, sigmaGrad);

    %define the single layer 'operator' object
    S=SingleLayer(kwave,Gamma);
    %construct Geometrical optics approximation on Gamma
    GOA=GeometricalOpticsApprox(uinc,Gamma);
    Sgoa=S*GOA;

    %get collocation points for each basis
    Xs = getColPoints( VHNAs );
    Xo = getColPoints( VHNAo );
    Xhp = getColPoints( Vhp );

    % initialise integration solvers
    %CG=CompGaussBasic(kwave,1000,1000);
    NSD=NSDlinearPhase(15,15);
    s=linspace(0,Gamma.L,100*kwave);
catch MatlabError
    display('An error occured in the problem setup');
    resultsData=strcat('Matlab error in setup: ',MatlabError.identifier);
    return;
end
% ------------------ solve for standard hp basis ------------------------
try
    preColMatrix=[];
    elCount=1;

    %create collocation matrix
    for b=Vhp.el
        Sb=S*b;
        preColMatrix=[preColMatrix Sb.col(Xhp)];
        elCount=elCount+1;
    end

    %create RHSa, no integrals here, so can be computed exactly
    RHS=f.eval(Xhp);

    % solve the integrals
    CG=CompGaussBasic(kwave,1000,1000);

    ColMatrix=zeros(size(preColMatrix));
    ColRHSb=zeros(length(Vhp.el),1);
    for n=1:length(Vhp.el)
        CG_=CG;
        fprintf('Computing standard column %d of %d\n',n,length(Vhp.el));
        for m=1:length(Vhp.el)
            ColMatrix(m,n)=CG_.eval(preColMatrix(m,n));
        end
    end
    %mix the RHSs, f(x)-A\Psi(x):
    coeffs_hp=ColMatrix\RHS;    
    v_hp=Projection(coeffs_hp,Vhp);
    clear preColMatrix coeffs_hp;
    
catch MatlabError
    display('An error occured in the standard hp mesh');
    resultsData=strcat('Matlab error in solving standard basis: ',MatlabError.identifier);
    return;
end
%---------------- solve for single mesh HNA space ------------------------
try 
    preColMatrix=[];
    elCount=1;
    tic;
    %create collocation matrix
    for b=VHNAs.el
        Sb=S*b;
        preColMatrix=[preColMatrix Sb.col(Xs)];
        elCount=elCount+1;
    end
    preColRHSb=Sgoa.col(Xs);
    ColMatrixNSD=zeros(size(preColMatrix));
    ColRHSb=zeros(length(VHNAs.el),1); ColRHSbNSD=ColRHSb;
    for n=1:length(VHNAs.el)
        fprintf('Computing SM column %d of %d\n',n,length(VHNAs.el));
        for m=1:length(VHNAs.el)
            ColMatrix(m,n)=NSD.eval(preColMatrix(m,n));
        end
        ColRHSbNSD(n)=NSD.eval(preColRHSb(n));
    end
    ts=toc;
    
    ColRHS=f.eval(Xs)-ColRHSbNSD;
    coeffs_s=ColMatrix\ColRHS;
    v_s=Projection(coeffs_s,VHNAo);
    err_s=sum(abs(v_s.eval(s)+GOA.eval(s)-v_hp.eval(s)))/sum(abs(v_hp.eval(s)));
    clear preColMatrix coeffs_s preColRHSb ColRHSbNSD;
    resultsData=char('CPU Time for single mesh HNA:',ts, 'Approx relative error:',err_s);
catch MatlabError
    display('An error occured in the single mesh');
    resultsData=strcat('Matlab error in solving single mesh HNA basis: ',MatlabError.identifier);
end%---------------- solve for overlapping mesh HNA space ------------------------
try 
    preColMatrix=[];
    elCount=1;
    tic;
    %create collocation matrix
    for b=VHNAs.el
        Sb=S*b;
        preColMatrix=[preColMatrix Sb.col(Xo)];
        elCount=elCount+1;
    end
    preColRHSb=Sgoa.col(Xo);
    ColMatrixNSD=zeros(size(preColMatrix));
    ColRHSb=zeros(length(VHNAo.el),1); ColRHSbNSD=ColRHSb;
    for n=1:length(VHNAo.el)
        fprintf('Computing OM column %d of %d\n',n,length(VHNAo.el));
        for m=1:length(VHNAo.el)
            ColMatrix(m,n)=NSD.eval(preColMatrix(m,n));
        end
        ColRHSbNSD(n)=NSD.eval(preColRHSb(n));
    end
    to=toc;
    clear preColMatrix;
    ColRHS=f.eval(Xs)-ColRHSbNSD;
    coeffs_o=ColMatrix\ColRHS;
    v_o=Projection(coeffs_o,VHNAo);
    err_o=sum(abs(v_o.eval(s)+GOA.eval(s)-v_hp.eval(s)))/sum(abs(v_hp.eval(s)));
    clear preColMatrix coeffs_o;
    resultsData=char('CPU Time for overlapping mesh HNA:',to, 'Approx relative error:',err_o);
catch MatlabError
    display('An error occured in the overlapping mesh');
    resultsData=strcat(resultsData,'Matlab error in solving single mesh HNA basis: ',MatlabError.identifier);
end
end

