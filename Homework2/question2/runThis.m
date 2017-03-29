%
%   CS 543 Homework 2 Question 2
%   Chongxin Luo
%   Feb. 19, 2017
%   Running files for Question 2

clc; clear all; close all;
obj1 = im2double(imread('\shape_align\object1.png'));
obj2 = im2double(imread('\shape_align\object2.png'));
obj2t = im2double(imread('\shape_align\object2t.png'));
obj3 = im2double(imread('\shape_align\object3.png'));

% Aline object2 to object2t
tic
[T1, outputimage1] = align_shape(obj2, obj2t);
toc
err1 = evalAlignment(outputimage1, obj2t);
dispim1 = displayAlignment(obj2, obj2t, outputimage1);
figure, colormap gray, imagesc(dispim1);
disp(sprintf('first aline error is : %s', err1 ));

% Aline object2 to object3
tic
[T2, outputimage2] = align_shape(obj2, obj3);
toc
err2 = evalAlignment(outputimage2, obj3);
dispim2 = displayAlignment(obj2, obj3, outputimage2);
figure, colormap gray, imagesc(dispim2);
disp(sprintf('second aline error is : %s', err2));

% Aline object2 to object1
tic
[T3, outputimage3] = align_shape(obj2, obj1);
toc
err3 = evalAlignment(outputimage3, obj1);
dispim3 = displayAlignment(obj2, obj1, outputimage3);
figure, colormap gray, imagesc(dispim3);
disp(sprintf('third aline error is : %s', err3));