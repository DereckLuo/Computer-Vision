function [T, outputimage] = align_shape(im1, im2)
% Function align shape: 
%   align image im1 and image m2, given a transformation T that maps
% non-zero points in im1 to non-zero points in im2
% Performing Iterative Closest Point (ICP) algorithm to find trandsform


% using matlab find to find all non-zero points to track in both image 
[Py, Px] = find(im1);
[Qy, Qx] = find(im2);
numPoints = size(Py, 1);
[imy, imx] = size(im2);

% perform initial transform of image im1 to line up with image im2
Pmeanx = mean(Px); Pmeany = mean(Py);
Pstdx = std(Px); Pstdy = std(Py);
Qmeanx = mean(Qx); Qmeany = mean(Qy);
Qstdx = std(Qx); Qstdy = std(Qy);

% figure, colormap gray, imagesc(im1);
% hold on 
% plot(Px(1), Py(1), 'r*');

transferMatrix = [1,0,Qmeanx;0,1,Qmeany;0,0,1]*([Qstdx,0,0; 0,Qstdy,0; 0,0,1]*([1/Pstdx,0,0; 0,1/Pstdy,0; 0,0,1]*([1,0,-Pmeanx; 0,1,-Pmeany; 0,0,1])));
T = transferMatrix;
for ii = 1:numPoints
    transpoint = transferMatrix * [Px(ii);Py(ii);1];
    Px(ii) = round(transpoint(1));
    Py(ii) = round(transpoint(2));
end

% initial transform check
% figure,colormap gray, imagesc(im2);
% hold on 
% plot(Px, Py, 'r.');
% Py(1)
% Px(1)

% variable used to track average distance as threshold 
oldDistance = 0;
newDistance = 100;

% repeat transformation prediction until reaching threshold 
while abs(oldDistance-newDistance) > 0.01
    % find nearest neighbor in image Q using bwdist()
    % D : Distance to nearst neightbor 
    % IDX : linear matrix index of nearest point
    [D, IDX] = bwdist(im2);

    % using Affine Transformations 
    % compute transformation T 
    X = []; y = [];
    distance = 0;
    for ii = 1:numPoints
        idx = IDX(Py(ii),Px(ii));
        distance = distance + D(Py(ii), Px(ii));
        qy = mod(idx, imy);
        qx = ceil(idx/imy);
        X = [X;
            Px(ii), Py(ii), 1, 0, 0, 0;
            0, 0, 0, Px(ii), Py(ii) 1];
        y = [y; double(qx); double(qy)];
    end 
    A = linsolve(X,y);
    Aff = [A(1), A(2), A(3);
           A(4), A(5), A(6);
           0   , 0   , 1];
    T = Aff * T;
    for ii = 1:numPoints 
        before = [Px(ii); Py(ii); 1];
        after = Aff * before;
        Px(ii) = round(after(1)); Py(ii) = round(after(2));
    end

    oldDistance = newDistance;
    newDistance = distance/numPoints;
end

% hold on 
% plot(Px, Py, 'g.');


outputimage = zeros(size(im1));
for ii = 1:numPoints
    outputimage(Py(ii), Px(ii)) = 1;
end 
%figure, colormap gray, imagesc(outputimage);
    
end 
