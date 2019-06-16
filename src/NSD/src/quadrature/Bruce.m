function [x, w] = Bruce( a,b,Npts,freq, singPoints, abWidth )
%brute force oscillatory integrator, used to validate NSD routines
p_max=12;
delta=0.15;

    if nargin<6
        abWidth=b-a;
    end
    
    if nargin<5
        singPoints=[];
    end
    
    if ~isempty(singPoints)
        copySingPoints=singPoints;
        clear singPoints;
        n=1;
        for j=1:length(copySingPoints)
            if a<=copySingPoints(j) && copySingPoints(j)<=b
                singPoints(n)=copySingPoints(j);
                n=n+1;
            end
        end
    end
    
    
            %first determine how many wavelengths are in interval
                wavelengths=ceil(abWidth/(2*pi/freq));
            x=[];
            w=[];
    
    if isempty(singPoints)
            %create some blank vectors
        %         x=zeros(wavelengths*Npts,1);
        %         w=x;
            %now fill them with composite Gaussian quadrature

            ChunkSize=abWidth/wavelengths;
                for j=1:wavelengths
                    [x_, w_] = quad_gauss(Npts, a+(j-1)*ChunkSize, a+j*ChunkSize);
                    x = [x; x_];
                    w = [w; w_];

                    %[x(((j-1)*Npts+1):(j*Npts)),w(((j-1)*Npts+1):(j*Npts))] = quad_gauss(Npts, a+(j-1)*abWidth/wavelengths, a+j*abWidth/wavelengths);
                end
                %now top it up to the end
        %         [x(((wavelengths-1)*Npts+1):end),w(((wavelengths-1)*Npts+1):end)] = quad_gauss(Npts, a+(wavelengths-1)*abWidth/wavelengths, b);
        return;
    else
        singNendPoints=unique([a singPoints b]);
        if length(singNendPoints)>2
            for s=1:(length(singNendPoints)-1)
                [x_, w_] = Bruce( singNendPoints(s),singNendPoints(s+1),Npts,freq, singPoints );
                x=[x; x_];
                w=[w; w_];
            end
        else
            [x_, w_]= GradedQuadFreq( Npts, p_max, delta,wavelengths );
            if isequal(singPoints,a)
                x=a+abWidth*x_;
                w=w_*abWidth;
            else %singPoints==b
                x=b-flipud(x_)*abWidth;
                w=(w_)*abWidth;
            end
        end
    end
end

