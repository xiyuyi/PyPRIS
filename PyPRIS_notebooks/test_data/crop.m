% crop
f1 = 'All_DH_NFz_NFx_NFy.tif'
f2 = 'Sample7_DH_FOV(00020).tif';
PSF = zeros(315, 3213, 180);
for i0 = 1 : 180
    im = imread(f1,'Index',i0);
    imwrite(im(4*63+[1:63],63*21+[1:63]), 'psf.tif', 'writemode', 'append');
end
for i0 = 1 : 31
    im = imread(f2,'Index',i0);
    imwrite(im(665:735,1275:1345), 'obsStack.tif', 'writemode', 'append');
end

