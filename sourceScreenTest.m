clc;
clear classes;

%wavenumber
kwave=128;

%define the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%inident plane wave
uinc=pointSource(kwave,[.5 1]);
    
%make an HNA basis on Gamma
pMax=5; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;
VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=singleLayer(kwave,Gamma);

[v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.15, 'progress');

%plot the output
s=linspace(0,Gamma.L,min(1000,10*kwave)).';
semilogy(s,abs(v_N.eval(s,1))./kwave); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([1E-7 1E4]); 
xlabel('s'); ylabel('|\partial u/\partial n - \partial u^i/\partial n|');

figure(2); plot(s,real(v_N.eval(s,1)+GOA.eval(s,1)),s,imag(v_N.eval(s,1)+GOA.eval(s,1))); 
xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([-.05*kwave .05*kwave]);
xlabel('s'); ylabel('\partial u/\partial n'); legend('real','imaginary');