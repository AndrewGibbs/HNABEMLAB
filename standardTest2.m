clc;
clear classes;
addPathsHNA();

%wavenumber
kwave=20;

%define the scatterer
vertices =   [0 0
              0 1];

% vertices =   [0     0;
%               1     0;];
      
%create 'edge' object for the screen
Gamma=edge(vertices);
% Gamma=edge(vertices);

% %inident plane wave
% %choose random angle:
% incAngle = rand*2*pi;
% uinc=planeWave(kwave,[cos(incAngle) sin(incAngle)]);

%inident plane wave
%uinc=planeWave(kwave,[-1 -1]/sqrt(2));

%uinc=planeWave(kwave,[0 -1]);
uinc=pointSource(kwave,[.5 1]);

%make an HNA basis on Gamma
pMax=4; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0; hMax=1/(2*kwave);
Vstd=hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);
DOFs=length(Vstd.el);
%define the single layer 'operator' object
S=singleLayer(kwave,Gamma);

tic;
[v_N, GOA, ~,~, colMatrix, colRHS] = ColHNA(S, Vstd, uinc, Gamma,'progress','oversample',1);
toc

%plot the output
% s=linspace(0,Gamma.L(1),20*kwave);
%plot(s,v_N.eval(s,1)/kwave);ylim([-kwave kwave]);

%now see how well an HNA basis fits to this:


[sALL, onSide] = getColPoints( Vstd, 5,1,'C');
for n = 1:Gamma.numSides
    %s=linspace(0,Gamma.L(n),min(1000,100*kwave)).';
    s = sALL(n==onSide);
    figure(n);
    fs{n} = v_N.eval(s,n)-GOA.eval(s,n);
    
    VHNA{n}=HNAsingleMesh(Gamma.side{n},pMax, kwave, throwAwayParam, nLayers, sigmaGrad,1);
    v_HNA{n} = VHNA{n}.leastSquares(s,fs{n});
    %plot(s,real(v_N.eval(s.',n)),s,imag(v_N.eval(s.',n)),s,real(GOA.eval(s.',n)),s,imag(GOA.eval(s.',n))); ylim([-2*kwave 2*kwave]);
    plot(s,real(fs{n}),'*',s,real(v_HNA{n}.eval(s,1)),'.'); ylim([-2*kwave 2*kwave]);
    legend('stdBEM','HNAbem');
end

HNAprojCoeffs = [v_HNA{1}.coeffs; v_HNA{2}.coeffs; v_HNA{3}.coeffs;];
stdCoeffs = v_N.coeffs;
save('stdProjHNAbasis','HNAprojCoeffs','stdCoeffs');
% s=linspace(0,Gamma.L,500*kwave);
% semilogy(s,abs(v_N.eval(s))./kwave); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([1E-4 1E7]); 
%VHNA.leastSquares(s,fs);