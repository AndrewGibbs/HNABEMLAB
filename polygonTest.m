clc;
clear classes;
addPathsHNA();

%wavenumber
kwave=100;%11*2*pi;

%define the screen
vertices =   [1    0;
              0     0;
              0    1];
      
%create 'edge' object for the screen/polygon
Gamma=polygon(vertices);

% %inident plane wave
% %choose random angle:
% incAngle = rand*2*pi;
% uinc=planeWave(kwave,[cos(incAngle) sin(incAngle)]);

%inident plane wave
uinc=planeWave(kwave,[-1 -1]./sqrt(2));
%uinc=planeWave(kwave,[0 -1]);

%make an HNA basis on Gamma
pMax=4; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;
VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=singleLayer(kwave,Gamma);

tic;
[v_N, GOA, colMatrix, colRHS, colMatrix2, colRHS2] = ColHNA(S, VHNA, uinc, Gamma,'oversample',1.2, 'progress');
toc

%plot the output
for n = 1:Gamma.numSides
    s=linspace(0,Gamma.L(n),min(1000,100*kwave));
    figure(n);
    plot(s,real(v_N.eval(s.',n)),s,imag(v_N.eval(s.',n))); ylim([-2*kwave 2*kwave]);
end

for n = 1:Gamma.numSides
    s=linspace(0,Gamma.L(n),20*kwave);
    figure(n);
    semilogy(s,abs(v_N.eval(s,n))./kwave); %ylim([1E-4 1E7]); 
end