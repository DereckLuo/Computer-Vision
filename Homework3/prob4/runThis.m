%%   CS 543 Homework 3 
%%   Chongxin Luo
%%   March. 19, 2017
%%   Running file for Question 4: Affine Structure from Motion

%% Part (a): form 3D structure 
clc; clear all; close all;
load('tracked_points.mat');
load('initial_keypoints.mat'); 

% m images n tracked points 
[m, n] = size(Xs);

% Step 1: for each frame, shift to centeriod of points for each camera (-
% mean)
for ii = 1:m
    meanx = mean(Xs(ii,:));
    for jj = 1:n
        Xs(ii,jj) = Xs(ii,jj) - meanx;
    end 
end
for ii = 1:m
    meany = mean(Ys(ii,:));
    for jj = 1:n
        Ys(ii,jj) = Ys(ii,jj) - meany;
    end 
end 
    
% form D matrix 
D = zeros(2*m, n);
idx = 1;
for ii = 1:m
    D(idx,:) = Xs(ii,:);
    D(idx+1,:) = Ys(ii,:);
    idx = idx + 2;
end 

% compute first estimation 
[U, W, V] = svd(D);

U3 = U(:, 1:3);
V3 = V(:, 1:3);
V3 = V3';
W3 = W(1:3, 1:3);

A = U3*sqrt(W3); X = sqrt(W3)*V3;

% Solve for orthographic constraints 
Astack = zeros(m*3, 9); 
Lstack = zeros(m*3, 1);
idx = 1; idx2 = 1;
for ii = 1:m
    A11 = A(idx,:); A12 = A(idx+1,:);
    a11 = A11(1); a12 = A11(2); a13 = A11(3); 
    a21 = A12(1); a22 = A12(2); a23 = A12(3);
    
    Astack(idx2  ,:) = [a11*a11, a12*a11, a13*a11, a11*a12, a12*a12, a13*a12, a11*a13, a12*a13, a13*a13];
    Astack(idx2+1,:) = [a21*a21, a22*a21, a23*a21, a21*a22, a22*a22, a23*a22, a21*a23, a22*a23, a23*a23];
    Astack(idx2+2,:) = [a21*a11, a22*a11, a23*a11, a21*a12, a22*a12, a23*a12, a21*a13, a22*a13, a23*a13];
    
    Lstack(idx2,:) = 1;
    Lstack(idx2+1,:) = 1;
    Lstack(idx2+2,:) = 0;
    
    idx = idx + 2;
    idx2 = idx2 + 3;
end 

L = Astack \ Lstack;
L = reshape(L, [3 3])';
C = chol(L, 'lower');
A = A * C;
X = inv(C)*X;


% plot in 3D point 
plot3(X(1,:), X(2,:), X(3,:), 'k.');


%% Part (b): find path for the camera
clc; close all;
I = zeros(size(A,1)/2,3);
J = zeros(size(A,1)/2,3);
idx = 1;
for ii = 1:size(A,1)/2
    I(ii,:) = A(idx,:);
    J(ii,:) = A(idx+1,:);
    idx = idx + 2;
end 

K = cross(I,J);

K = normr(K);

figure, plot(K(:,1),1:m,'r*');
figure, plot(K(:,2),1:m,'g*');
figure, plot(K(:,3),1:m,'b*');
