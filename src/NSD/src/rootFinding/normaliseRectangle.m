function rectOut = normaliseRectangle( rect,  height, width)
%surprisingly difficult process of converting a complex rectangle into a
%normalised square
if nargin==1
     height=1;
     width=1;
end
    rectOutR=real(rect)./abs(real(rect));
    rectOutL= imag(rect)./abs(imag(rect)) ;
    rectOutR(isnan(rectOutR))=0;
    rectOutL(isnan(rectOutR))=0;
    rectOut=width*rectOutR+height*1i*rectOutL;

end

