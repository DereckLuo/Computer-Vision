clc; clear all; close all;
temp = ceil(7/3)
temp2 = mod(11853,344)
%%
% solving matrix
clc; clear all; close all;
ori = [2;3;1]
trans = [1,2,3;4,5,6;0,0,1]
dest = trans*ori
solve = linsolve(dest, ori)
%% Affine prediction test
clc; clear all; close all;
x = [1,2,1,0,0,0;
    0,0,0,1,2,1];
y = [5;4];
x
y
c = linsolve(x,y)
f = [c(1),c(2),c(3);
     c(4),c(5),c(6);
     0,0,1]
 ret = f * [1;2;1]
 %%
tic
x1 = [1,2,3;
      4,5,6;
      7,8,9]
x2 = [10,11,12;
      13,14,15;
      16,17,18]
x3 = x1*x2
y = [1;2;3]
ret1 = x2*y;
ret1 = x1*ret1
ret2 = x3*y
toc