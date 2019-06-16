function [yn,P, nearCP] = nearlyFiniteCheck(p0, h0, h1, startCP, CPs, freq)
%a condition to determine it a path is nearly finite

    jumpThresh=50;
    distThresh1=1/freq;
    distThresh2=.25; % don't ever get this close to another critical point
    
    p=p0(2:end);
   
    N=length(h1);
    otherCPs = CPs(startCP~=CPs);
    diff_h1=(abs(abs(h1(1:(N-1)))-abs(h1(2:N)))) ./ (p(2:N)-p(1:(N-1)));
    if ~isempty(otherCPs)
        dist_h1=abs(otherCPs-h0);
    else
        dist_h1=inf;
    end
    
    
    yn=false;
    P=[];
    nearCP=[];
        
    %first condition checks for jumps in derivative, second checks that
    %path doesn't get too close to stationary point. It can start close
    %though.
%     if (max(diff_h1)>jumpThresh && min(min(dist_h1))<distThresh1) || min(min(dist_h1))<distThresh2 && min(min(dist_h1))~=min(dist_h1(1,:))
%         yn=true;
        if min(min(dist_h1))<distThresh2 && min(min(dist_h1))~=min(dist_h1(1,:))
            yn=true;
            dist=inf;
            for CPindex=1:length(CPs)
                if ~ismember(startCP,CPs(CPindex))
                    [dist_, pIndex_] = min(abs(CPs(CPindex)-h0));
                    if dist_<dist
                        nearCP=CPindex;
                        dist=dist_;
                        pIndex=pIndex_;
                    end
                end
            end
            P=p(pIndex);
        elseif max(diff_h1(2:end))>jumpThresh && min(min(dist_h1))<distThresh1
            yn=true;
            [~,pIndex] = max(diff_h1);
            P=p(pIndex);
            dist=inf;
            for CPindex=1:length(CPs)
                %only check other critical points
                if ~ismember(startCP,CPs(CPindex))
                    dist_=abs(h0(pIndex)-CPs(CPindex));
                    if dist_<dist
                        nearCP=CPindex;
                        dist=dist_;
                    end
                end
            end
        end
        %return;
%     else
%     end
    
end