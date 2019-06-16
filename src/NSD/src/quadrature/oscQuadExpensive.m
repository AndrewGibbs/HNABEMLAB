function [x,w] = oscQuadExpensive(a,b,freq, Npts, abWidth)
%brute force quadrature for path without singularities
            %first determine how many wavelengths are in interval
            if nargin==4
                abWidth=abs(b-a);
            end
            abDir=(b-a)/abWidth;
            wavelengths=ceil(abWidth/(2*pi/freq));
            x=[];
            w=[];
            ChunkSize=abWidth/wavelengths;
            [x_, w_] = quad_gauss(Npts, 0, ChunkSize);
            for j=1:wavelengths
                %[x_, w_] = quad_gauss(Npts, a+(j-1)*ChunkSize, a+j*ChunkSize);
                x = [x; a+(j-1)*ChunkSize*abDir + x_*abDir];
                w = [w; w_*abDir];

                %[x(((j-1)*Npts+1):(j*Npts)),w(((j-1)*Npts+1):(j*Npts))] = quad_gauss(Npts, a+(j-1)*abWidth/wavelengths, a+j*abWidth/wavelengths);
            end
end

