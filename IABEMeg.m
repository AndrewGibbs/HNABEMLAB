clc;
clear classes;

%wavenumbers
K = 2.^(13:18);
std_kmax =1030;

%define the screen
vertices =   [0    0;
              1    0];
          
sourcePoint = [.5 1]     ;     
      
%create 'screen' object
Gamma=edge(vertices);
    
%make an HNA basis on Gamma
pMax=4; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;

n=0;
tNSD = NaN(size(K));
tCG = NaN(size(K));
for kwave = K
    fprintf('\nk = %d', kwave);
    n = n+1;
    %inident plane wave
    %uinc=planeWave(kwave,[1 1]./sqrt(2));
    uinc=pointSource(kwave,sourcePoint);%planeWave(kwave,sqrt(2)*[-1 -1]);
    
    %make HNA basis
    VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
    DOFs=length(VHNA.el);
    %define the single layer 'operator' object
    S=singleLayer(kwave,Gamma);

    %solve using pathFinder quadrature
    fprintf('\nPathFinder quadrature');
    tic;
    [v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.25, 'progress');
    tNSD(n) = toc
    if kwave<std_kmax
        %solve using standard BEM
        fprintf('\nStandard quadrature');
        tic;
        [v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.25, 'progress','CG');
        tCG(n) = toc
    end
    
save('IABEMdata_','K','tNSD','tCG');
end
loglog(K,tNSD,'x',K,tCG,'x');
xlabel('wavenumber');
ylabel('CPU time (seconds)');
legend('Steepest descent quadrature','Composite Gauss quadrature');