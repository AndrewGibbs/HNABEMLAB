function [h,dhdp, W, FPindices] = getBranchesOfFinitePaths(prePathLengthsVec, pathEndIndex, criticalPoints, pathPowers, G, freq, N, COVsingularities, RelTol)

 
%FPindices=false(length(prePathLengthsVec),2);

%the condition which has flagged up the finite paths at this point is
%necessary, but not sufficient. Some of these paths will not be finite. In
%such a case, the endpoints of the two components of the finite path will
%be really far apart, and would be obvious (to a human). If the two components of the path
%are farther apart than the following constant, we assume they are not
%connected.
maxJoinDistThresh = .1;%.1;

    %initialise output
    for j=1:length(prePathLengthsVec)
        for m=1:pathPowers(j)
            h{j,m}=[];
            dhdp{j,m}=[];
            W{j,m}=[];
            FPindices{j}(m)=0;
        end
    end
    for j=1:length(prePathLengthsVec)
       if ~isinf(prePathLengthsVec(j))
          %there's a finite path here
          P1=(prePathLengthsVec(j)/2)^(1/pathPowers(j));
          ICsSD = NSDpathICv3( pathPowers(j), G, criticalPoints(j), false);
          [Psd{j}, Wsd{j}] = pathQuadFinite( P1, COVsingularities, pathPowers(j), freq, N );
          for m=1:pathPowers(j)
               [~,Hsd{m}] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(j)-1,G, ICsSD{m}, false), [0; Psd{j}], ICsSD{m}, odeset('RelTol',RelTol) ); 
              %avoid the zero/IC entry, as this didn't correspond to a
              %quadrature point
              Hsd{m}=Hsd{m}(2:end,:);     
          end
          
          P2=(prePathLengthsVec(j)-P1^(pathPowers(j)))^(1/pathPowers(pathEndIndex(j)));
          ICsSA = NSDpathICv3( pathPowers(pathEndIndex(j)), G, criticalPoints(pathEndIndex(j)), true);
          [Psa{j}, Wsa{j}] = pathQuadFinite( P2, COVsingularities, pathPowers(pathEndIndex(j)), freq, N );
          for n=1:pathPowers(pathEndIndex(j))
              [~,Hsa{n}] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(pathEndIndex(j))-1,G, ICsSA{n}, true), [0; Psa{j}], ICsSA{n}, odeset('RelTol',RelTol) );
              %avoid the zero/IC entry, as this didn't correspond to a
              %quadrature point
              Hsa{n}=Hsa{n}(2:end,:);
              %flip this so endpoints are at the end of joined path
              Hsa{n}=flipud(Hsa{n});
          end
       
       %now check which was the best match. See which midpoints of which
       %branches end up the closest
       dist=inf;
           for  m=1:pathPowers(j)
               for n=1:pathPowers(pathEndIndex(j))
                    dist_=abs(Hsd{m}(end,1)-Hsa{n}(1,1));
                    if dist_<dist
                        dist=dist_;
                        mClosest=m;
                        nClosest=n;
                    end
               end
           end
           
           if dist<maxJoinDistThresh
                h{j,mClosest}=[Hsd{mClosest}(:,1); Hsa{nClosest}(:,1);];
                dhdp{j,mClosest}=[Hsd{mClosest}(:,2); Hsa{nClosest}(:,2);];
                W{j,mClosest}=[Wsd{j}; Wsa{j};];
                FPindices{j}(mClosest)=pathEndIndex(j);
           end
           clear ICsSA Hsd Hsa Wsa Wsd dist;
       end
    end
end