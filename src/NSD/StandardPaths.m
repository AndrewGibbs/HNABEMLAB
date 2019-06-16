function StandardPaths()

%adds all NSD45 relevant paths:
%addpath NSD;
% addpath NSD/pathChoosing;
% addpath NSD/pathFinding;
% addpath NSD/quadrature;
% addpath NSD/rootFinding;
% addpath NSD/singularPoints;
% addpath NSD/DaveCode


addpath src;
addpath src/pathChoosing;
addpath src/pathFinding;
addpath src/quadrature;
addpath src/rootFinding;
addpath src/singularPoints;
addpath src/infContours;

%some of Daan's old hoipack paths get used too:
% addpath hoipack;
% addpath hoipack/expansion;
% addpath hoipack/quadrules;
% addpath hoipack/OPQ;

%Moved Daan's hoipack stuff into quadrature folder:
addpath src/quadrature/expansion;
addpath src/quadrature/quadrules;
addpath src/quadrature/OPQ;

%(should probably go through these at some point and delete what's not
%needed)

end