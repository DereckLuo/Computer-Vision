%loading image
clc;
close all;
clear all;
levels = 5;
sigma = 2;
hsize = 3;
I = im2double(imread('tower.jpg'));
I = rgb2gray(I);
G_filter = fspecial('gaussian', hsize, sigma);

G_pyramid = cell([1,levels]);
L_pyramid = cell([1,levels]);
G_fft = cell([1,levels]);
L_fft = cell([1,levels]);
%imshow(I);
%Matlab buildin Gaussian pyramid function for testing
%for i = 1:levels
%    I = impyramid(I, 'reduce');
%    figure, colormap winter, imagesc(log(abs(fftshift(fft2(I)))));
%    figure, imshow(I);
%end

%---- Generating both G pyramid and L pyramid ------%
for n = 1:levels
    prev_img = I; 
    %post_img = imresize(imgaussfilt(I,sigma),0.5);
    post_img = imresize(imfilter(I,G_filter), 0.5);
    %post_img = impyramid(I, 'reduce');
    up_scale = imresize(post_img, size(prev_img));
    Laplace = prev_img - imgaussfilt(up_scale, sigma);
    G_fft{n} = log(abs(fftshift(fft2(I)))+1);
    L_fft{n} = log(abs(fftshift(fft2(Laplace)))+1);
    G_pyramid{n} = I;
    L_pyramid{n} = Laplace;
    I = post_img;
end
%---- output both pyramid ----%
figure 
ha = tight_subplot(2,levels,[.02 .02],[.02 .02], [.05, .05]);
for n = 1:levels; axes(ha(n)); colormap gray; imagesc(G_pyramid{n}); end
for n = 1:levels; axes(ha(n+levels)); colormap gray; imagesc(L_pyramid{n}); end
set(ha(2:1),'XTickLabel',''); set(ha(2:1),'YTickLabel','')

figure 
ta = tight_subplot(2,levels,[.02 .02],[.02 .02], [.05, .05]);
for n = 1:levels; axes(ta(n)); colormap jet; imagesc(G_fft{n}, [0,10]); end
for n = 1:levels; axes(ta(n+levels)); colormap jet; imagesc(L_fft{n}, [0,10]); end
set(ta(2:1),'XTickLabel','');set(ta(2:1),'YTickLabel','')




