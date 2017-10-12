P=[1:8];
K=[64 128 256 512 1024 2048 4096];
SampleScale=linspace(1,2,10);
%now plot stuff.
k_=0;
lines={'b.-','ko-','rx-','g+-','c*-','ms-','yd-','bv-'};
LEG={'p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8'};

%hold on;
for k=K
    k_=k_+1;
    %load(sprintf('k%dp8_hpBEM',k));
    figure(k_);
    p_=0;
    for p=P
        p_=p_+1;
        %semilogy(P,err(P,k_,s_));
        %do this next bit manually because matlab can't tell a vector from
        %a 
        s_=0;
        for s=SampleScale
            s_=s_+1;
            yData(s_)=err(p_,k_,s_);
            CPsPerDOF(s_)=CPs(p_,k_,s_)/DOFs(p_,k_,s_);
        end
        %semilogy(SampleScale,yData,lines{p_});
        semilogy(CPsPerDOF,yData,lines{p_});
        title(sprintf('k=%d',k));
        hold on;
    end
    legend(LEG)
    xlabel('ColPoints / DOFs')
    ylabel('Relative error')
    %ylim([10^-2.5,10^.5]);
    hold off;
end