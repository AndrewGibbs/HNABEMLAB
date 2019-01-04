# Preliminaries & version info

Fast solver for high frequency scattering problems in 2-dimensions. Uses collocation BEM with an Hybrid Numerical Asymptotic (HNA) method, which absorbs oscillatory part of solution into approximation space.

Requires PathFinder (available from https://github.com/AndrewGibbs/NSDpackage) to be on Matlab search path.

For polygons (not yet thoroughly tested), requires Chebfun (https://github.com/chebfun/chebfun) to be on Matlab search path.

Currently seems stable for screens with plane wave incidence. Can (in principle) also hand polygons, and point source incidence. Computes all of these in frequency independent time, using an HNA basis [1] and oscillatory quadrature routines [3], with oversampled collocation [2].

# Problem statement & formulation

For problems of scattering of an incident wave <img src="http://latex.codecogs.com/svg.latex?u^i=\mathrm{e}^{\mathrm{i}k\mathbf{d}\cdot\mathbf{x}}" border="0"/> by a screen <img src="http://latex.codecogs.com/svg.latex?\Gamma" border="0"/>, we aim to compute the total field <img src="http://latex.codecogs.com/svg.latex?u:=u^i+u^s" border="0"/>, where <img src="http://latex.codecogs.com/svg.latex?u^s" border="0"/> is the scattered field.

Solution satisfies the Helmholtz equation:

![equation](https://latex.codecogs.com/gif.latex?%28%5CDelta&plus;k%5E2%29u%3D0%5Cquad%5Ctext%7Bin%20%7D%5Cmathbb%7BR%7D%5E2%5Csetminus%5CGamma%2C%20%5Cquad%20u%3D0%5Cquad%5Ctext%7Bon%20%7D%5CGamma)

which we can reformulate for <img src="http://latex.codecogs.com/svg.latex?[\partial_nu]:=\partial_n^+u - \partial_n^-u" border="0"/> as

![equation](https://latex.codecogs.com/gif.latex?%5Cint_%5CGamma%5Cfrac%7B%5Cmathrm%7Bi%7D%7D%7B4%7DH_0%5E%7B%281%29%7D%28k%7C%5Cmathbf%7Bx%7D-%5Cmathbf%7By%7D%7C%29%5B%5Cpartial_nu%5D%28%5Cmathbf%7By%7D%29%5Cmathrm%7Bd%7Ds%28%5Cmathbf%7By%7D%29%3Du%5Ei%28%5Cmathbf%7Bx%7D%29%2C%5Cquad%5Ctext%7Bon%20%7D%5CGamma.)

The HNA Ansatz for the screen is

 <img src="http://latex.codecogs.com/svg.latex?[\partial_nu](\mathbf{x})(s):=v_+(s)\mathrm{e}^{\mathrm{i}ks}+v_-(s)\mathrm{e}^{-\mathrm{i}ks}+\Psi(\mathbf{x}),\quad\text{on }\Gamma" border="0"/>

where ![equation](https://latex.codecogs.com/gif.latex?%5CPsi%3D%5Cfrac%7B%5Cpartial%20u%5Ei%7D%7B%5Cpartial%20n%7D), and <img src="http://latex.codecogs.com/svg.latex?v_\pm" border="0"/> are non-oscillatory [1].

# Approxixmation

We approximate <img src="http://latex.codecogs.com/svg.latex?=v_+(s)\mathrm{e}^{\mathrm{i}ks}+v_-(s)\mathrm{e}^{-\mathrm{i}ks}" border="0"/> using an HNA basis on a single mesh on <img src="http://latex.codecogs.com/svg.latex?\Gamma" border="0"/>, graded towards the endpoints to capture the singularities.

We solve our discrete system using an oversampled collocation method, as outlined in [2], taking around 40% more collocation points than basis elements. This leads to a rectangular system, which can be solved in a least-squares sense, via a truncated SVD.

# Example code
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

And plot the solution in the domain (this bit isn't frequency independent):

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

# References:

[1] D Hewett, S. Langdon & S.N. Chandler-Wilde, <a href="https://arxiv.org/pdf/1401.2786.pdf">A frequency-independent boundary element method for scattering by two-dimensional screens and apertures</a>

[2] B. Adcock & D. Huybrechs, <a href="https://arxiv.org/pdf/1802.01950.pdf">Frames and numerical approximation II: generalized sampling</a>

[3] A. Gibbs & D. Huybrechs, <a href="https://people.cs.kuleuven.be/~andrew.gibbs/AGibbs5.pdf">Frames and numerical approximation II: generalized sampling</a>
