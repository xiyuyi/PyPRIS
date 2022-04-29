% trim  PSF
f1 = 'psf_BP+250.tif';
f2 = 'psf_BP-250.tif';
c1 = 91;
c2 = 41;

f = f1;
c = c1;
for i0 = 1:151
    im = imread(f, 'Index', i0);
    imwrite(im-100,['rmbg_',f],'writemode','append');
end


f = f2;
c = c2;
for i0 = 1:151
    im = imread(f, 'Index', i0);
    imwrite(im-100,['rmbg_',f],'writemode','append');
end

