# Preliminaries & version info

Fast solver for high frequency scattering problems in 2-dimensions. Uses collocation BEM with an Hybrid Numerical Asymptotic (HNA) method, which absorbs oscillatory part of solution into approximation space.

Currently stable for screens with plane wave incidence, and multiple alligned screens, which are not too close together. Can (in principle) also handle polygons, and point source incidence, although this has not been analysed and tested, and is not fully frequency independent.

# Problem statement & formulation

This is outline *very* briefly here, for a full explanation, please see [1].

For problems of scattering of an incident wave
$$u^i(\mathbf{x})=\mathrm{e}^{\mathrm{i}k\mathbf{d}\cdot\mathbf{x}}$$
 by a screen $\Omega=\Gamma$, or polygon $\Omega$ with boundary $\Gamma$, we aim to compute the total field $u:=u^i+u^s$, where $u^s$ is the scattered field.

Our solution satisfies the exterior Helmholtz BVP:

$$
(\Delta+k^2)u=0\quad\text{in}\quad\mathbb{R}^2\setminus\overline\Omega
$$
 
 and $u^s$ satisfies the Sommerfeld radiation condition.

For simplicity, we focus on the screen problem in this explanation. This can be reformulated [1] as $[\partial_nu]:=\partial_n^+u - \partial_n^-u$ as

$$\frac{\mathrm{i}}{4}\int_\Gamma^{~} H_0^{(1)}(k|\mathbf{x}-\mathbf{y}|)[\partial_nu](\mathbf{y})\mathrm{d}s(\mathbf{y}) = u^i(\mathbf{x}),\quad\text{on}\quad\Gamma$$

The HNA Ansatz for the screen is

$$[\partial_nu](\mathbf{x})(s):=v_+(s)\mathrm{e}^{\mathrm{i}ks}+v_-(s)\mathrm{e}^{-\mathrm{i}ks}+\Psi(\mathbf{x}),\quad\text{on }\quad\Gamma"
$$

where  $\Psi:=2\partial_n^+u^i$, and $v_\pm$ are non-oscillatory [1].

HNABEMLAB uses this ansatz (and a similar approach for the screen), and works by approximating $v_\pm$.

# Numerical Approxixmation

First create instances of the fundamental objects which define our problem.

```matlab
    % run addPathsHNA() to add necessary search paths
    %wavenumber
    kwave = 60;

    %create 'screen' object ---------------------------------------------------
    vertices =   [0    0;
                  1    0];
    Gamma = Screen(vertices);

    %inident plane wave -------------------------------------------------------
    d = [1 1]./sqrt(2); %direction as a vector
    uinc = planeWave(kwave, d);
```
As is standard for HNA methods, we approximate 

$$\nu_N(s)\approx~[\partial_nu](\mathbf{x}(s))-\Psi(\mathbf{x}(s))=v_+(s)\mathrm{e}^{\mathrm{i}ks}+v_-(s)\mathrm{e}^{-\mathrm{i}ks}$$

using an HNA basis on a single mesh on $\Gamma$ graded towards the endpoints to capture the singularities.  Hence our solution

$$\nu_N\in\{\rho_+(s)\mathrm{e}^{\mathrm{i}ks}+\rho_-(1-s)\mathrm{e}^{-\mathrm{i}ks}:\rho_\pm\in\mathbb{P}_{p,n}(0,1)\}$$

where $\mathbb{P}_{p,n}(0,1)$ is the space of piecewise polynomials on the unit interval of order $p$ with $n$ layers of grading towards zero.

```matlab

    %make an HNA basis on Gamma -----------------------------------------------
    pMax = 8; %polynomial degree
    cL = 2; %layers of grading per polynomial degree
    sigmaGrad = 0.15; %grading ratio
    nLayers = cL*(pMax+1)-1; %number of layers of grading
    throwAwayParam = 0; %no need to remove any basis elements
    OverSample = 1.4; %choose amount to oversample by (40% here)
    % construct the HNA basis (single mesh):
    VHNA = HNAsingleMesh(Gamma, pMax, kwave, throwAwayParam, nLayers, sigmaGrad, 1);
    DOFs = length(VHNA.el); %get total #DOFs

    % construct the single layer potential 'operator' ---------------------------
    S=singleLayer(kwave,Gamma);
```

We solve our discrete system using an oversampled collocation method, as outlined in [2], taking around 40% more collocation points $\{s_m\}_{m=1}^M$ than basis elements. This leads to a rectangular system, which can be solved in a least-squares sense, via a truncated SVD, minimising

$$\sum_{m=1}^M|\mathcal{S}_k(\nu_N(s_m)+\Psi(\mathbf{x}(s_m)))-u^i(\mathbf{x}(s_m))|\quad\text{for}\quad~M\geq1.4N$$

```matlab
    %solve (and time)
    tic;
    [v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample', OverSample, 'progress');
    T = toc
    
    
```
Now plot the solution on the boundary:

```matlab

    s = linspace(0,Gamma.L,min(1000,10000*kwave));
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

The code is also frequency independent for multiple screens, by instantiating the object

```matlab
  Gamma=MultiScreen(vertices,[.0 .2 .4 .6 .8 1]);
```

where the vector in the second argument corresponds to points at which the two vertices are split. Below is the same example for multiple alligned screens, with $k=100$

![HNABEMLAB](https://raw.github.com/AndrewGibbs/HNABEMLAB/master/mutliScreenPlot.png)

# References:

[1] A. Gibbs, D Hewett, D. Huybrechs, E. Parolin <a href="https://doi.org/10.1007/s42985-020-00013-3">Fast hybrid numerical-asymptotic boundary element methods for high frequency screen and aperture problems based on least-squares collocation, SN Partial Differ. Equ. Appl, 1 (2020): 1-26.
