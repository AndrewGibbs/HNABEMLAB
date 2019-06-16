clc;
clear classes;

%have just tried to simplify & generalise code, so that polygons, screens
%and multiple screens all fit inside of the same framework.

%define some different geometries
sqV = [1 1; 1 0; 0 0; 0 1;];
scV = [0 0; 1 0;];
scS = [0 0.3 0.5 1];

singleScreen = Screen(scV);
twoScreens = MultiScreen(scV, scS);
square = ConvexPolygon(sqV);

%define an incident wave
kwave = 50;
d = [1 -1]./sqrt(2);
ui = planeWave(kwave,d);

%now define all types of basis on each type of scatterer
pMax = 4;
hMax = 1/kwave;
nLayers = pMax;
sigmaGrad = .15;
scats = {singleScreen, twoScreens, square};
for n=1:3
    tensorBasis{n,1} = HNAoverlappingMesh(scats{n}, pMax, kwave, nLayers, sigmaGrad);
    tensorBasis{n,2} = HNAsingleMesh(scats{n}, pMax, kwave, 0, nLayers, sigmaGrad);
    tensorBasis{n,3} = hpStandardBasis(scats{n}, pMax, hMax, nLayers, sigmaGrad);
end

%now define GOA on each type
for n=1:3
    GOA{n} = GeometricalOpticsFunction(ui,scats{n});
end

V = tensorBasis{2,3};

Sk=singleLayer(kwave,V.obstacle);
ColHNAv2(Sk, V, ui, V.obstacle);