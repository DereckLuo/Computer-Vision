% Function getKeypoints: 
%  given a input image, get the keypoints for tracking using second moment
%  matrix 
%
function [keyXs, keyYs] = getKeypoints(im, tau)

scale = 0.03; tau = 0.002
%im = im2double(imread('\tracking\images\hotel.seq0.png'));
%im = im2double(imread('Capture.JPG'));
%im = rgb2gray(im);
% smooth the image and remove noise for gradient calculation
im = imgaussfilt(im);

% compute gradient in x and y direction
[Gx, Gy] = imgradientxy(im); 

% compute square of derivatives
Ix2 = Gx.^2;
Iy2 = Gy.^2;
IxIy = Gx.*Gy;

% Gaussian filter again to sum up the local neighborhood 
GIx2 = imgaussfilt(Ix2, 2);
GIy2 = imgaussfilt(Iy2, 2);
GIxIy = imgaussfilt(IxIy, 2);

% compute harris score 
har = (GIx2.*GIy2) - (GIxIy.^2) - scale.*((GIx2+GIy2).^2);

% perform non-max suppression on 5x5 window with threshold tau
[keyXs, keyYs] = nonmax(har, tau);

end

