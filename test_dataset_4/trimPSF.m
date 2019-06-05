% trim PSF stack
f1 = 'psf_BP+250.tif'
f2 = 'psf_BP-250.tif'
f1sta = 59;
f1end = 136;
f2sta = 14;
f2end = 80;



for i0 = 1:f1sta-1
    imwrite(zeros(61,61,'uint16'),['trimmed_',f1],'writemode','append');
end
for i0 = f1sta:f1end
    im = imread(f1, 'Index', i0);
    imwrite(im,['trimmed_',f1],'writemode','append');
end
for i0 = f1end-1:151
    imwrite(zeros(61,61,'uint16'),['trimmed_',f1],'writemode','append');
end

for i0 = 1:f2sta-1
    imwrite(zeros(61,61,'uint16'),['trimmed_',f2],'writemode','append');
end
for i0 = f2sta:f2end
    im = imread(f2, 'Index', i0);
    imwrite(im,['trimmed_',f2],'writemode','append');
end
for i0 = f2end-1:151
    imwrite(zeros(61,61,'uint16'),['trimmed_',f2],'writemode','append');
end