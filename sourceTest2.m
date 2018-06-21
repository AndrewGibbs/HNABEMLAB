clc;
clear classes;

%wavenumber
kwave=128;

%define the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%inident point source
uinc=pointSource(kwave,[.5 1]);
    
%make an HNA basis on Gamma
pMax=4; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;
VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=singleLayer(kwave,Gamma);

tic;
[v_CG, GOA, colMatrix_CG, colRHS_CG] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.15, 'progress', 'CG');
toc

tic;
[v_NSD, GOA, colMatrix_NSD, colRHS_NSD] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.15, 'progress' );
toc

matErr = (abs(colMatrix_CG - colMatrix_NSD)./abs(colMatrix_CG));
RHSerr = (abs(colRHS_CG - colRHS_NSD)./abs(colRHS_CG));

%plot the output
s=linspace(0,Gamma.L,500*kwave);
figure(1); semilogy(s,abs(v_CG.eval(s,1))./kwave); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([1E-4 1E7]); 
figure(2); plot(s,real(v_CG.eval(s,1))); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([-2*kwave 2*kwave]); 