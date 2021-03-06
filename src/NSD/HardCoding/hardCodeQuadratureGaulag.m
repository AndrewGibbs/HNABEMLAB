
name = 'gausLagHC';
FID = fopen(strcat(name,'.m'),'w');
N = 100;

fprintf(FID,'function [x,w] = %s(n)\n',name);

fprintf(FID,'switch n\n');

for n=1:N
    fprintf(FID,'case %d\n',n);
    [x,w] = quad_gauss_laguerre(n);
    fprintf(FID,'\tx = [');
    for m = 1:length(x)
        fprintf(FID,'%.16f ',x(m));
    end
    fprintf(FID,"].';\n");
    fprintf(FID,'\tw = [');
    for m = 1:length(x)
        fprintf(FID,'%.16f ',w(m));
    end
    fprintf(FID,"].';\n");
end
fprintf(FID,'otherwise\n');
fprintf(FID,'[x,w] = quad_gauss_laguerre(n);\n');
fprintf(FID,'end');