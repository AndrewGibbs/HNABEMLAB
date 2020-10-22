clc;
clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
kwave=20;

%create 'screen' object ---------------------------------------------------
vertices =   [0    0;
              1    0];
Gamma=Screen(vertices);

%inident plane wave -------------------------------------------------------
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 4; %polynomial degree
cL = 2; %layers of grading per polynomial degree
hMax = 2*pi/(2*kwave);
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any b    asis elements
OverSample = 1.25; %choose amount to oversample by (40% here)
% construct the HNA basis (single mesh):
%VHNA = HNAsingleMesh(Gamma, pMax, kwave, throwAwayParam, nLayers, sigmaGrad, 1);
%VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
Vstd = hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);

% construct the single layer potential 'operator' ---------------------------
S=singleLayer(kwave,Gamma);

%solve
[v_h, GOA] = ColHNA(S, Vstd, uinc, Gamma,'oversample', OverSample, 'progress', 'symmetry');

%plot the output
s=linspace(0,Gamma.L,min(10000,20*kwave));
s = s(2:(end-1));

n=1; %one side on screen
s=linspace(0, Gamma.component(n).L, 20*kwave);
plot(s,real(v_h.eval(s, n)),'k', s, real(GOA.edgeComponent(n).eval(s.')),'r');
%plot(s,real(GOA.edgeComponent(n).eval(s.')));
ylim([-2*kwave 2*kwave]);
legend('v_h', 'GOA');
return;

figure(1);
%semilogy(s,fliplr(abs(v_N.eval(s,1))),s,fliplr(abs(v_N.evalComp(s,'+',1))),':',s,fliplr(abs(v_N.evalComp(s,'-',1))),'--','LineWidth',2);
plot(s,fliplr((v_N.eval(s,1))),s,abs(fliplr((v_N.evalComp(s,'+',1)))),'r:',s,abs(fliplr((v_N.evalComp(s,'-',1)))),'g--','LineWidth',2); ylim([-20 20]);
xlabel('s');
hold on;
plot(s,-abs(fliplr((v_N.evalComp(s,'+',1)))),'r:',s,-abs(fliplr((v_N.evalComp(s,'-',1)))),'g--','LineWidth',2); ylim([-20 20]);
leg = legend('Re\{v_N(s)\}','\pm|v^+_{1,N}(s)|','\pm|v^-_{1,N}(s)|');
xlim([Gamma.component.supp(1)-.1 Gamma.component.supp(2)+.1] );
set(gca,'FontSize',16);
leg.FontSize = 16;
title(sprintf('k=%d',kwave),'FontSize',16);
%legend('Re(v_N(s))','Re(v^+_{1,N}(s))','Re(v^-_{1,N}(s))');
figure(2);
plot(s,mod(angle(fliplr((v_N.evalComp(s,'+',1)))),2*pi),':',s,mod(angle(fliplr((v_N.evalComp(s,'-',1)))),2*pi),'--','LineWidth',2);
leg = legend('arg(v^+_{1,N}(s))','arg(v^-_{1,N}(s))');
xlabel('s');
xlim([Gamma.component.supp(1)-.1 Gamma.component.supp(2)+.1] );
set(gca,'FontSize',16);
leg.FontSize = 16;
title(sprintf('k=%d',kwave),'FontSize',16);
%return;
%plot the far-field pattern:
figure(3);
theta = linspace(0,2*pi,50*kwave);
Fv_N = FarField_lessSlow_stillSteady(Gamma, v_N, kwave, theta);
FPsi = FarField_lessSlow_stillSteady(Gamma, GOA, kwave, theta);
semilogy(theta,abs(FPsi+Fv_N));
return;%'fontsize',

%now plot the solution in the domain:
figure(3);
domainPlot(Gamma,uinc,GOA,v_N,kwave);
axis equal;

