function [P, endIndex, finitePathPowers, A] = finitePathTest( SPs, G, pathPowers, thresh )
%function takes in stationary points as 'SPs' and returns matrix of path
%lengths, if finite

    if nargin==3
        threshBase=1E-6;
    end
    N=length(SPs);
    P=inf(N,1);
    finitePathStart=false(N,1);
    finitePathEnd=false(N,1);
    endIndex=NaN(N,1);
    finitePathPowers=NaN(N,2);
    g=G{1};

    for n=1:N
       for m=1:N
          %scale threshold by distance between stationary points
          thresh = threshBase * abs(SPs(n)-SPs(m));
          test=(g(SPs(m))-g(SPs(n)))/1i;
          if SPs(m)==SPs(n)
              A(m,n)=0; %the following condition contains some serious bodging.
              %let's talk about what each bit does:
              %1. checks the real part is constaint along finite path
              %2&3 checks the argument of the path h(p) is positive real
          elseif abs(real(g(SPs(m)))-real(g(SPs(n))))<thresh && real(test)>-thresh && abs(imag(test))<thresh
              A(m,n)=real(test);
          else
              A(m,n)=inf;
          end
       end
    end
    
    %in case it's possible that one finite path crosses two stationary
    %points, the following will truncate the path length at the first:
    for n=1:N
        for m=1:N
           if A(m,n)<P(n) && 0<A(m,n)
                P(n)=A(m,n);
                finitePathStart(n)=true;
                finitePathEnd(m)=true;
                endIndex(n)=(m);
                finitePathPowers(n,1)=pathPowers(n);
                finitePathPowers(n,2)=pathPowers(m);
           end
        end
    end

end