%
%   CS 543 Homework 2 Question 1
%   Chongxin Luo
%   Feb. 19, 2017
%   Running files for Question 1

clc; clear all; close all;
I = im2double(imread('\tracking\images\hotel.seq0.png'));

% part (1.1) getting the key points in the image
% setting parameter threshold tau = 0.002
tau = 0.002;
[keyXs, keyYs] = getKeypoints(I, tau);
% display keypoints and red *
figure, colormap gray
imagesc(I)
hold on
plot(keyXs, keyYs, 'r*');

% part (1.2) tracking all keypoints over 50 frams 
% display random 20 points tracking, point out points moved out of frame 
tracking(keyXs, keyYs);


% %%
% % using matlab given corner function to double check result
% %clc; clear all; close all;
% I = im2double(imread('\tracking\images\hotel.seq0.png'));
% colormap gray
% C = corner(I);
% figure, color cray
% imagesc(I)
% hold on 
% plot(C(:,1),C(:,2),'r*');
