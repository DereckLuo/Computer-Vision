wavelength = 20;
orientation = [0 45 90 135];
g = gabor(wavelength,orientation);
imshow(real(g(4).SpatialKernel),[]);
%%
%%
clc;
close all;
clear all;
orientations = [0:30:180];
fColumn = length(orientations);
G_filter = cell([1,fColumn]);
for ii = 1:fColumn
    G_filter{ii} = customgauss([21,21],5,2,orientations(ii),0,1,[0,0]);
    add = sum(sum(G_filter{ii}));
    G_filter{ii} = G_filter{ii}/add;
end
ha = tight_subplot(2,4,[.02 .02],[.02 .02], [.05 .05]);
for n = 1:7; axes(ha(n)); imagesc(G_filter{n}); end
%%
clc;
close all;
clear all;
ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

%%
clc; close all; clear all;
I = im2double(imread('sample2.jpg'));
I = rgb2gray(I);
figure 
imshow(I);
figure
g = gabor(2, 90);
output = imgaborfilt(I,g);
imshow(output);
%%
clc; close all; clear all;
I = im2double(imread('sample2.jpg'));
I = rgb2gray(I);
G_filter = fspecial('prewitt');
G_filter2 = rot90(G_filter);
output = imfilter(I,G_filter);

imagesc(output);

figure
imagesc(imfilter(I,G_filter2));

figure
G_filter3 = [1,0,-1.5;
             0 1,0;
             -1.5,0,1]
imagesc(imfilter(I,G_filter3));
%%
clc; close all; clear all;
G_filter = fspecial('gaussian', 4*3+1, 3);
colorbar
imagesc(G_filter);
colorbar
G_filter
sum(sum(G_filter))
I = im2double(imread('sample2.jpg'));
figure
imagesc(I);
figure
imagesc(imfilter(I,G_filter));
%%
clc; close all; clear all;
sigx = 3;
sigy = 1;
width = 13;
orientations = [0:30:180];
fColumn = length(orientations);

G_filter = cell([1, fColumn]);


for ii = 1:fColumn
    G_filter{ii} = customgauss([width,width],sigx, sigy,orientations(ii),0,1,[0,0]);
    add = sum(sum(G_filter{ii}));
    G_filter{ii} = G_filter{ii}/add;
    figure
    imagesc(G_filter{ii});
end
%%
clc; close all; clear all;
I = im2double(imread('sample2.jpg'));
F = fspecial('gaussian', 13, 3)
width = 13;
sigx = 3;
sigy = 1;
G_filter = customgauss([width,width],sigx, sigy,45,0,1,[0,0]);
imagesc(G_filter)
%%
clc; close all; clear all;
I = im2double(imread('sample.jpg'));
G = fspecial('gaussian',6,2);
I = imfilter(I,G);
figure
imagesc(I);
figure
F = fspecial('laplacian', 0.2);
output = imfilter(I, F);
imagesc(output)
