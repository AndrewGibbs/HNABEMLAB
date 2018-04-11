clc;
clear classes;
messageFlag=1;
%wavenumber
kwave=128;

%define the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%inident plane wave
uinc=planeWave(kwave,[1 1]./sqrt(2));
    
%make an HNA basis on Gamma
pMax=2; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;
VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=singleLayer(kwave,Gamma);

tic;
[v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.15);
toc

%plot the output
s=linspace(0,Gamma.L,500*kwave);
semilogy(s,abs(v_N.eval(s))./kwave); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([1E-4 1E7]); 