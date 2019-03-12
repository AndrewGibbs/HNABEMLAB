function [ X, t, onSide, W] = getColPoints( Vbasis, overSamplesPerMeshEl, scaler, type, side )
%returns collocation points that have been evenly spread across basis
%elements
 %ChebyshevRoots( 3, 'Tn', [1 2] )
 
 if nargin == 4
     side = 1;
 end
 
 if isa(Vbasis.obstacle,'polygon')
%      %call self recursively
%      t = [];
%      tR = []; %new variable, distance to far corner of edge, to avoid rounding errors
%      onSide = [];
%      side = [];
     X = [];
     onSide = [];
     t =[];
     for n = 1:Vbasis.obstacle.numSides
         [X_,t_,~] = getColPoints( Vbasis.edgeBasis{n}, overSamplesPerMeshEl, scaler, type, n );
%          t = [t; t_];
%          tma = [tma; tma_];
%          bmt = [bmt; bmt_];
         onSide((length(onSide)+1):(length(onSide)+length(t_))) = n;
         X = [X X_];
         %onSide = [onSide; onSide_];
         t = [t; t_];
     end
     return;
 end
 
 
 if nargin <=1
     overSamplesPerMeshEl=1;
 end
 if nargin<=2
     scaler=1;
 end
 if overSamplesPerMeshEl<1
     error('Oversampling param must be at least one');
 end
 
 if nargin <=3
    type=('C');%chebyshev
 end
 
    t = []; tma = []; bmt = [];
    if isa(Vbasis,'HNAoverlappingMesh')
%         E=[Vbasis.mesh{1}.el Vbasis.mesh{2}.el];
%         M=[Vbasis.meshDOFs{1} Vbasis.meshDOFs{2}];
          [E_, M] = Vbasis.mimicSingleMesh;
          E = E_.el;
          overlapElFlag = true;
    else
        E=Vbasis.mesh.el;
        M=Vbasis.meshDOFs;
        midEl = [];
        overlapElFlag = false;
    end
    sCount = 0;
    for m=1:length(M)
        %pts=ceil((m.pMax+1)*overSamplesPerMeshEl);
        pts=ceil(M(m)*overSamplesPerMeshEl);
        if strcmp(type,'C')%Chebyshev
            s=sort(cos(pi*(2*(1:pts)-1)/(2*pts))).';
            %w =ones(pts,1)*pi./pts;
            w = RiemannWeights(s,-1,1);
        elseif strcmp(type,'U')%uniform
            s=linspace(-1,1,pts+2).'; %add two extra points (endpoints)
            s=s(2:(end-1)); %delete the extra two points
            w = ones(pts,1)/pts;
        else
            error('Point distribution type not recognised, must be uniform or Chebyshev');
        end
        
        %now scale everything onto the mesh element:
        t  = [t; E(m).interval(1)+E(m).width*.5*(s+1);];
        %w  = [w; E(m).width*w_*.5];
        %tR = [tR; E(m).distR + E(m).width*.5*(flipud(s)+1);];
        tma = [tma; E(m).width*.5*(s+1);];
        bmt = [bmt; E(m).width*.5*(flipud(s)+1);];
        sSubCount = 0;
        for s_=s.'
            sCount = sCount + 1;
            sSubCount = sSubCount+1;
%             if m == midEl
%                 overlapElFlag = true;
%             else
%                 overlapElFlag = false;
%             end
            X(sCount) = collocationPoint( E(m), .5*(s_+1), side, m, E(m).width*w(sSubCount)*.5, overlapElFlag);
        end
    end
    onSide = ones(length(t),1);
end