function X = getColPointsV2( Vbasis, overSamplesPerMeshEl, type )
%returns collocation points that have been evenly spread across basis
%elements

    xCount = 0;
    for component=1:Vbasis.obstacle.numComponents
         if nargin <=1
             overSamplesPerMeshEl=1;
         end
         if overSamplesPerMeshEl<1
             error('Oversampling param must be at least one');
         end

         if nargin <=3
            type=('C');%chebyshev
         end

        t = []; tma = []; bmt = [];
        if isa(Vbasis,'HNAoverlappingMesh')
              [E_, M] = Vbasis.mimicSingleMesh;
              E = E_.el;
              overlapElFlag = true;
        else
            E=Vbasis.mesh{component}.el;
            M=Vbasis.edgeBasis{component}.meshDOFs;
            midEl = [];
            overlapElFlag = false;
        end
        
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
            tma = [tma; E(m).width*.5*(s+1);];
            bmt = [bmt; E(m).width*.5*(flipud(s)+1);];
            xSubCount = 0;
            for s_=s.'
                xCount = xCount + 1;
                xSubCount = xSubCount+1;
                X(xCount) = collocationPoint( E(m), .5*(s_+1), component, m, E(m).width*w(xSubCount)*.5, overlapElFlag);
            end
        end
    end
end