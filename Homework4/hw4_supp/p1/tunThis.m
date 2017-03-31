clc; close all; clear all;
im = im2double(imread('house2.jpg'));
[cIndMap, time, imgVis] = slic(im, 100, 10);