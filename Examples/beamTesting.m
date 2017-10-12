%changes to make:
%-%decide on a convention for domain ordering, and adjust code accordingly
%-%should put a check in SingleLayerR2 to make sure function is on right
%domain, similar with L2 inner products

clear classes;
kwave=1;

%define the screen
vertices=[0 0       %first vertex
          2*pi 0];     %second vertex
      
%create 'edge' object for the screen
Gamma1=edge(vertices);
Gamma2=edge(vertices+3); %shift it a bit for the second obstacle

%define a density on the second obstacle
f=density(Gamma1,@(s) 1i);
g=density(Gamma2,@(s) 1);

Sg=SingleLayerR2(kwave,g,Gamma2);

Gg=GeometricalOpticsApprox(Sg,Gamma1);

%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma1);

%example of operator S acting on operator G acting on function g
SGg=S*Gg;

%example of 3D inner product
I=L2(SGg,f);
CG=CompGaussBasic(kwave,1000,20); 
L2_of_SGg_and_g=CG.eval(I)
% J=L2(Gg,f);
% L2_of_Gg_and_g=CG.eval(J);