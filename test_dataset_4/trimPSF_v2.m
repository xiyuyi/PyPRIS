% trim  PSF
f1 = 'psf_BP+250.tif';
f2 = 'psf_BP-250.tif';
c1 = 91;
c2 = 41;

f = f1;
c = c1;
for i0 = 1:c-41
    imwrite(ones(61,61,'uint16').*100,['trimmed2_',f],'writemode','append');
end
for i0 = c-40:c+40
    im = imread(f, 'Index', i0);
    imwrite(im,['trimmed2_',f],'writemode','append');
end
for i0 = c+41:151
    imwrite(ones(61,61,'uint16').*100,['trimmed2_',f],'writemode','append');
end


f = f2;
c = c2;
for i0 = 1:c-41
    imwrite(zeros(61,61,'uint16'),['trimmed2_',f],'writemode','append');
end
for i0 = c-40:c+40
    im = imread(f, 'Index', i0);
    imwrite(im,['trimmed2_',f],'writemode','append');
end
for i0 = c+41:151
    imwrite(ones(61,61,'uint16').*100,['trimmed2_',f],'writemode','append');
end
