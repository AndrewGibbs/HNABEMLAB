function t = getColPoints( Vbasis, overSamplesPerMeshEl, scaler, type )
%returns collocation points that have been evenly spread across basis
%elements
 %ChebyshevRoots( 3, 'Tn', [1 2] )
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
 
    t=[];
    if isa(Vbasis,'HNAoverlappingMesh')
        E=[Vbasis.mesh{1}.el Vbasis.mesh{2}.el];
        M=[Vbasis.meshDOFs{1} Vbasis.meshDOFs{2}];
    else
        E=Vbasis.mesh.el;
        M=Vbasis.meshDOFs;
    end
    for m=1:length(M)
        %pts=ceil((m.pMax+1)*overSamplesPerMeshEl);
        pts=ceil(M(m)*overSamplesPerMeshEl);
        if strcmp(type,'C')%Chebyshev
            s=sort(cos(pi*(2*(1:pts)-1)/(2*pts))).';
        elseif strcmp(type,'U')%uniform
            s=linspace(-1,1,pts+2).'; %add two extra points (endpoints)
            s=s(2:(end-1)); %delete the extra two points
        else
            error('Point distribution type not recognised, must be uniform or Chebyshev');
        end
        if mod(length(s),2)==0 %if even number of points
            mid=length(s)/2;
            s1=s(1:mid);
            s2=[];
            s3=s((mid+1):end);
        else %odd number
            mid=(length(s)+1)/2;
            s1=s(1:(mid-1));
            s2=s(mid);
            s3=s((mid+1):end);
        end
        %now scale everything towards endpoints
        s1=-1+(s1+1).*scaler;
        s3=1-(1-s3).*scaler;
        s=[s1;s2;s3;];
        %now scale everything onto the mesh element:
        t=[t; E(m).interval(1)+E(m).width*.5*(s+1);];
    end
end

