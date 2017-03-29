%%   CS 543 Homework 3 
%%   Chongxin Luo
%%   March. 19, 2017
%%   Running file for Question 3: Epipolar Geometry

%% Part (a1): Solving for Fundemental Matrix F
clc; clear all; close all;
im1 = im2double(imread('chapel00.png'));
im2 = im2double(imread('chapel01.png'));
load('prob3.mat');

% normalize image coordinates 
mean1x = mean(c1); mean1y = mean(r1);
mean2x = mean(c2); mean2y = mean(r2);

std1 = 0;
for ii = 1:size(c1,1)
    std1 = std1 + sqrt((c1(ii)-mean1x)*(c1(ii)-mean1x) + (r1(ii)-mean1y)*(r1(ii)-mean1y));
end
std1 = std1/size(c1,1);
std2 = 0;
for ii = 1:size(c2,1)
    std2 = std2 + sqrt((c2(ii)-mean2x)*(c2(ii)-mean2x) + (r2(ii)-mean2y)*(r2(ii)-mean2y));
end
std2 = std2/size(c2,1);
% std1x = std(c1); std1y = std(r1);
% std2x = std(c2); std2y = std(r2);

% forming T matrix (homogenious x,y,1)
% Tp = [1/std1x, 0, -mean1x/std1x; 
%       0, 1/std1y, -mean1y/std1y;
%       0, 0, 1];
T1 = [1/std1, 0, -mean1x/std1;
      0, 1/std1, -mean1y/std1;
      0, 0, 1];
T2 = [1/std2, 0, -mean2x/std2;
      0, 1/std2, -mean2y/std2;
      0, 0, 1];
% % Tq = [1/std2x, 0, -mean2x/std2x;
% %       0, 1/std2y, -mean2y/std2y;
% %       0, 0, 1];

% normalize points and represent into homogenious coordinates 
norm1 = [c1, r1, ones(size(c1,1),1)];
norm2 = [c2, r2, ones(size(c2,1),1)];
for ii = 1:size(c1,1)
    normpoint = T1*[norm1(ii,1);norm1(ii,2);norm1(ii,3)];
    norm1(ii,1) = normpoint(1);
    norm1(ii,2) = normpoint(2);
    norm1(ii,3) = normpoint(3);
end
for ii = 1:size(c2,1)
    normpoint = T2*[norm2(ii,1);norm2(ii,2);norm2(ii,3)];
    norm2(ii,1) = normpoint(1);
    norm2(ii,2) = normpoint(2);
    norm2(ii,3) = normpoint(3);
end


totalpoints = size(matches,1);
threshold = 1.1;
bestF = []; best_inliers = []; best_outliers = [];
targetRatio = 0.5;
bestRatio = 0;
while bestRatio < targetRatio
    % initialize containers 
    inliers = []; outliers = [];
    
    % random select 8 points 
    pointset = randperm(totalpoints, 8);
    
    % Building matrix F
    A = zeros(8,9);
    for ii = 1:8
        p = matches(pointset(ii),1); q = matches(pointset(ii),2);
        u = norm1(p,1); v = norm1(p,2);
        up = norm2(q,1); vp = norm2(q,2);
        A(ii, :) = [u*up, u*vp, u, v*up, v*vp, v, up, vp, 1];
    end
    
    % solve f using SVD 
    [U,S,V] = svd(A);
    f = V(:, end);
    F = reshape(f, [3 3])';
    % resolve det(F) = 0 to constraint using SVD
    [U,S,V] = svd(F);
    S(3,3) = 0;
    F = U*S*V';
    % De-normalize F
    normF = transpose(T1)*F*T2;
    
    % Test F on all points 
    for ii = 1:totalpoints
        point1 = [c1(matches(ii,1)); r1(matches(ii,1)); 1];
        point2 = [c2(matches(ii,2)); r2(matches(ii,2)); 1];
        
        l1 = normF * point2;
        distance1 = abs(l1(1)*point1(1) + l1(2)*point1(2) + l1(3))/sqrt(l1(1)*l1(1)+l1(2)*l1(2));
        
        l2 = normF' * point1;
        distance2 = abs(l2(1)*point2(1) + l2(2)*point2(2) + l2(3))/sqrt(l2(1)*l2(1)+l2(2)*l2(2));
        
        if distance1 <= threshold || distance2 <= threshold
            inliers = [inliers; ii];            
        else 
            outliers = [outliers; ii];
        end 
    end 
    
    % restore best result so far 
    bestRatio = size(inliers,1)/totalpoints;
    if bestRatio >= targetRatio
        bestF = normF;
        best_inliers = inliers;
        best_outliers = outliers;
    end
end

disp(bestF);
bestnormF = bestF/norm(bestF);
disp(norm(bestF, 2));
disp(norm(bestnormF));

% plot outliers 
figure, colormap gray, imagesc(im1)
hold on
for ii = 1:size(best_outliers,1)
    idx = matches(best_outliers(ii),1);
    plot(c1(idx), r1(idx), 'g.');
end

% plot seven random points and their epipolar lines 
pointset = randperm(size(best_inliers,1), 7);
figure, colormap gray, imagesc(im1);
hold on 
for ii = 1:7
    idx = best_inliers(ii);
    p1 = [c1(matches(idx,1)); r1(matches(idx,1)); 1];
    p2 = [c2(matches(idx,2)); r2(matches(idx,2)); 1];
    l1 = bestF*p2;
    
    plot(p1(1), p1(2), 'r+');
    leftx = 0; rightx = 500;
    lefty = (l1(1)*leftx + l1(3))/(-l1(2));
    righty = (l1(1)*rightx + l1(3))/(-l1(2));
    plot([leftx, rightx], [lefty, righty], 'g');
end

figure, colormap gray, imagesc(im2);
hold on
for ii = 1:7
    idx = best_inliers(ii);
    p1 = [c1(matches(idx,1)); r1(matches(idx,1)); 1];
    p2 = [c2(matches(idx,2)); r2(matches(idx,2)); 1];
    l2 = bestF'*p1;
    
    plot(p2(1), p2(2), 'r+');
    leftx = 0; rightx = 500;
    lefty = (l2(1)*leftx + l2(3))/(-l2(2));
    righty = (l2(1)*rightx + l2(3))/(-l2(2));
    plot([leftx, rightx], [lefty, righty], 'g');
end



  




