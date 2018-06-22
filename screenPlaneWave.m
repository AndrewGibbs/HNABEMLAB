clc;
clear classes;

%wavenumber
kwave=11*2*pi;

%define the screen
vertices =   [0    1;
              1    0];
      
%create 'screen' object
Gamma=edge(vertices);

%inident plane wave
%uinc=planeWave(kwave,[1 1]./sqrt(2));
uinc=planeWave(kwave,[0 -1]);
    
%make an HNA basis on Gamma
pMax=4; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;
VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=singleLayer(kwave,Gamma);

tic;
[v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.25, 'progress');
toc

%plot the output
for n = 1:Gamma.numSides
    s=linspace(0,Gamma.L(n),min(1000,10*kwave));
    figure(n);
    semilogy(s,abs(v_N.eval(s,n))./kwave); ylim([1E-4 1E7]); 
end
% s=linspace(0,Gamma.L,500*kwave);
% semilogy(s,abs(v_N.eval(s))./kwave); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([1E-4 1E7]); 