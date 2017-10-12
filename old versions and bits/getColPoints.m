function s = getColPoints( Vbasis, overSamplesPerMeshEl, scaler )
%returns collocation points that have been evenly spread across basis
%elements
 %ChebyshevRoots( 3, 'Tn', [1 2] )
 if nargin <=1
     overSamplesPerMeshEl=0;
 end
%  if nargin <=2
%         AlrPar = @(a,b,abWidth,N) (a+abWidth*.5*(1+cos(pi*(2*(1:N)-1)/(2*N)))).';
%  else
%      AlrPar = @(a,b,abWidth,N) (a+abWidth*.5*(1+max(1,(1-scaler)/(cos(pi/(2*N))))*cos(pi*(2*(1:N)-1)/(2*N)))).';
%  end
 %scaler2 = @(N) 2*(1-scaler)/(1+cos(pi/(2*N)));
    AlrPar = @(a,b,abWidth,N) (a+abWidth*.5*(1+cos(pi*(2*(1:N)-1)/(2*N)))).';
 
    s=[];
    if isa(Vbasis,'HNAoverlappingMesh')
        for meshEl = Vbasis.mesh{1}.el
            %s=[s; ChebyshevRoots( meshEl.pMax+1+overSamplesPerMeshEl, 'Tn', meshEl.interval ).';];
            s_=AlrPar(meshEl.interval(1),meshEl.interval(2),meshEl.width,meshEl.pMax+1+overSamplesPerMeshEl);
            if mod(length(s),2)==0 %if even number of points
            s=[s; AlrPar(meshEl.interval(1),meshEl.interval(2),meshEl.width,meshEl.pMax+1+overSamplesPerMeshEl);];
        end
        for meshEl = Vbasis.mesh{2}.el
            %s=[s; ChebyshevRoots( meshEl.pMax+1+overSamplesPerMeshEl, 'Tn', meshEl.interval ).';];
            s=[s; AlrPar(meshEl.interval(1),meshEl.interval(2),meshEl.width,meshEl.pMax+1+overSamplesPerMeshEl);];
        end
    elseif isa(Vbasis,'hpStandardBasis')
        for meshEl = Vbasis.mesh.el
            s=[s; ChebyshevRoots( meshEl.pMax+1+overSamplesPerMeshEl, 'Tn', meshEl.interval ).';];
        end
    elseif isa(Vbasis,'HNAsingleMesh')
        for meshEl = Vbasis.mesh.el
            s=[s; ChebyshevRoots( (meshEl.pMax+1+overSamplesPerMeshEl)*(meshEl.osc+1), 'Tn', meshEl.interval ).';];
        end
    end
end

