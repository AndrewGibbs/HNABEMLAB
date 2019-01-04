# Preliminaries & version info

Fast solver for high frequency scattering problems in 2-dimensions. Uses collocation BEM with an Hybrid Numerical Asymptotic (HNA) method, which absorbs oscillatory part of solution into approximation space.

Requires PathFinder (available from https://github.com/AndrewGibbs/NSDpackage) to be on Matlab search path.

For polygons, requires Chebfun (https://github.com/chebfun/chebfun) to be on Matlab search path.

Currently stable for screens with plane wave incidence. Can (in principle) also hand polygons, and point source incidence. Computes all of these in frequency independent time, using an HNA basis and oscillatory quadrature routines, with oversampled collocation.

# Problem statement & formulation

For problems of scattering of an incident wave <img src="http://latex.codecogs.com/svg.latex?u^i=\mathrm{e}^{\mathrm{i}k\mathbf{d}\cdot\mathbf{x}}" border="0"/> by a screen <img src="http://latex.codecogs.com/svg.latex?\Gamma" border="0"/>, we aim to compute the total field <img src="http://latex.codecogs.com/svg.latex?u:=u^i+u^s" border="0"/>, where <img src="http://latex.codecogs.com/svg.latex?u^s" border="0"/> is the scattered field.

Solution satisfies the Helmholtz equation:

![equation](https://latex.codecogs.com/gif.latex?%28%5CDelta&plus;k%5E2%29u%3D0%5Cquad%5Ctext%7Bin%20%7D%5Cmathbb%7BR%7D%5E2%5Csetminus%5CGamma%2C%20%5Cquad%20u%3D0%5Cquad%5Ctext%7Bon%20%7D%5CGamma)

which we can reformulate as

![equation](https://latex.codecogs.com/gif.latex?%5Cint_%5CGamma%5CPhi%28x%2Cy%29%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20n%7D%28y%29dy%3Du%5Ei%28x%29%2C%5Cquad%5Ctext%7Bon%20%7D%5CGamma.)
```matlab
    % run addPathsHNA() to add necessary search paths
    %wavenumber
    kwave=60;

    %create 'screen' object ---------------------------------------------------
    vertices =   [0    0;
                  1    0];
    Gamma=edge(vertices);

    %inident plane wave -------------------------------------------------------
    d = [1 1]./sqrt(2); %direction as a vector
    uinc=planeWave(kwave,d);

    %make an HNA basis on Gamma -----------------------------------------------
    pMax = 8; %polynomial degree
    cL = 2; %layers of grading per polynomial degree
    sigmaGrad=0.15; %grading ratio
    nLayers = cL*(pMax+1)-1; %number of layers of grading
    throwAwayParam = 0; %no need to remove any basis elements
    OverSample = 1.4; %choose amount to oversample by (40% here)
    % construct the HNA basis (single mesh):
    VHNA = HNAsingleMesh(Gamma, pMax, kwave, throwAwayParam, nLayers, sigmaGrad, 1);
    DOFs = length(VHNA.el); %get total #DOFs

    % construct the single layer potential 'operator' ---------------------------
    S=singleLayer(kwave,Gamma);

    %solve (and time)
    tic;
    [v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample', OverSample, 'progress');
    T = toc
    
    
```
Now plot the solution on the boundary:

```matlab

    s=linspace(0,Gamma.L,min(1000,10000*kwave));
    figure;
    semilogy(s,abs(v_N.eval(s,1))); ylim([1E-2 1E3]);
    xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] );
    ylabel('|\partial u/\partial n - \Psi|')

```
![HNABEMLAB](https://raw.github.com/AndrewGibbs/HNABEMLAB/master/boundaryPlot_k60.png)

And plot the solution in the domain:

```matlab
    figure;
    domainPlotPoints = min(1000,kwave*20);
    y = linspace(-1,1,domainPlotPoints);
    x = linspace(-1,2,domainPlotPoints);
    Sv = singleLayerDomain(Gamma, v_N, kwave, x, y);
    SPsi = singleLayerDomain(Gamma, GOA, kwave, x, y);
    [X1, X2] = meshgrid(x,y);
    u_N = uinc.eval(X1,X2) - Sv - SPsi;
    imagesc(x,y,real(u_N))
    shading interp;
    hold on;
    Gamma.draw;
```
![HNABEMLAB](https://raw.github.com/AndrewGibbs/HNABEMLAB/master/domainPlot_k60.png)
