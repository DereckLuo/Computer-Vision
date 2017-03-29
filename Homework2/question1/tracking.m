function tracking(keyXs, keyYs)
%Function tracking: given a set of tracking points
% read 50 frame images, and track the movement

dir = '\tracking\images\';
pointSize = size(keyXs);
numPoints = pointSize(2);
%selectXs = []; selectYs = [];
selectXs = keyXs; selectYs = keyYs;
storeXs = selectXs'; storeYs = selectYs'; outarr = [];
load('tracking\initial_keypoints.mat');
% random select 20 points 
randnum = randi(numPoints,1,20);

I = im2double(imread('\tracking\images\hotel.seq0.png'));
figure, colormap gray, imagesc(I)
hold on
plot(selectXs, selectYs, 'r*');
for ii = 0:49
    disp(ii);
    file0 = sprintf('hotel.seq%d.png',ii);
    file1 = sprintf('hotel.seq%d.png',ii+1);
    path0 = fullfile(dir, file0);
    path1 = fullfile(dir, file1);

    im0 = im2double(imread(path0));
    im1 = im2double(imread(path1));
    [selectXs, selectYs, out] = predictTranslationAll(selectXs, selectYs, im0, im1);
    outarr = [outarr, out];
    storeXs = [storeXs, selectXs'];
    storeYs = [storeYs, selectYs'];
end

hold on
% plot(storeXs, storeYs, 'g.', 'linewidth', 3);
for ii = 1:size(randnum,2)
    idx = randnum(ii);
    plot((storeXs(idx,:)), (storeYs(idx,:)), 'g.', 'linewidth', 3);
end

for ii = 1:size(outarr,2)
    idx = outarr(ii);
    plot(keyXs(idx), keyYs(idx), 'b*', 'linewidth', 3);
end

end

