function [cIndMap, time, imgVis] = slic(img, K, compactness)

%% Implementation of Simple Linear Iterative Clustering (SLIC)
%
% Input:
%   - img: input color image
%   - K:   number of clusters
%   - compactness: the weighting for compactness
% Output: 
%   - cIndMap: a map of type uint16 storing the cluster memberships
%   - time:    the time required for the computation
%   - imgVis:  the input image overlaid with the segmentation
clc; clear all; close all;
img = im2double(imread('house2.jpg'));
K = 100; compactness = 10;

tic;
% Put your SLIC implementation here

% initial variables 
[y, x, z] = size(img); % image x and y size 
S = ceil(sqrt((x*y)/K));  % grid space S 
disp(S);
% clusters = zeros(K, 5); % K x 5 array to store clusters 
clusters = [];
disp(size(img)); 
% initial cluster centers
idx = 1;
for ii = 1:S:y
    for jj = 1:S:x
        if idx > K
            break;
        end
        % if ii <= y && jj <= x
        % x_pos = int16(jj); y_pos = int16(ii);
        % disp([idx, x_pos, y_pos]);
        % temp = [jj, ii, img(ii,jj,1), img(ii,jj,2), img(ii,jj,3)];
        clusters = [clusters; [jj, ii, img(ii,jj,1), img(ii,jj,2), img(ii,jj,3)]];%clusters : x, y, r, g, b
        % clusters(idx, :) = [jj, ii, img(ii,jj,1), img(ii,jj,2), img(ii,jj,3)]; 
        idx = idx + 1;  
        %end
    end 
end 
numpix = zeros(size(clusters,1),1);

% print initial cluster points 
imagesc(img);
hold on
for ii = 1:size(clusters,1)
    plot(clusters(ii,1), clusters(ii,2), 'r.');
end

% move cluster points to smallest gradient
[Gmag, Gdir] = imgradient(rgb2gray(img));
for ii = 1:size(clusters,1)
    xcord = clusters(ii,1); ycord = clusters(ii,2);
    minGrad = Gmag(ycord, xcord);
    minx = xcord; miny = ycord;
    disp([xcord, ycord]);
    for n = max(1,ycord-1) : min(ycord+1,y)
        for k = max(1,xcord-1) : min(xcord+1,x)
            curGrad = Gmag(n,k);
            if curGrad < minGrad
                minGrad = curGrad; minx = k; miny = n;
            end 
        end 
    end 
  clusters(ii,1) = minx; clusters(ii,2) = miny;
end 

% print initial cluster points 
figure, imagesc(img);
hold on
for ii = 1:size(clusters,1)
    plot(clusters(ii,1), clusters(ii,2), 'r.');
end

% initialize pixel labels and pixel distances
labels = zeros(y,x); distances = zeros(y,x);
labels(:) = -1; distances(:) = inf;

% repeat to compute distance and move cluster centers
for l = 1:50
    disp(norm(clusters));
    for ii = 1:size(clusters,1)
        xcord = clusters(ii,1); ycord = clusters(ii,2);
        for n = max(1,round(ycord-S)) : min(round(ycord+S),y)
            for k = max(1,round(xcord-S)): min(round(xcord+S),x)
                r1 = img(n,k,1); g1 = img(n,k,2); b1 = img(n,k,3);
                r2 = clusters(ii,3); g2 = clusters(ii,4); b2 = clusters(ii,5);
                dc = sqrt((r1-r2)*(r1-r2)+(g1-g2)*(g1-g2)+(b1-b2)*(b1-b2));
                ds = sqrt((xcord-k)*(xcord-k)+(ycord-n)*(ycord-n));
                D = sqrt(dc*dc + (ds/S)*(ds/S)*compactness*compactness);
                if D < distances(n,k)
                    distances(n,k) = D;
                    labels(n,k) = ii;
                end
            end
        end
    end

    % recompute cluster center
    clusters = zeros(size(clusters));
    numpix = zeros(size(clusters,1),1);
    for ii = 1:y
        for jj = 1:x
            clusternum = labels(ii,jj);
            r = img(ii,jj,1); g = img(ii,jj,2); b = img(ii,jj,3);
            clusters(clusternum,1) = clusters(clusternum,1) + jj;
            clusters(clusternum,2) = clusters(clusternum,2) + ii;
            clusters(clusternum,3) = clusters(clusternum,3) + r;
            clusters(clusternum,4) = clusters(clusternum,4) + g;
            clusters(clusternum,5) = clusters(clusternum,5) + b;
            numpix(clusternum) = numpix(clusternum) + 1;
            % clusters(clusternum,6) = clusters(clusternum,6) + 1;
        end 
    end 
    for ii = 1:size(clusters,1)
%         clusters(ii,1) = clusters(ii,1)/clusters(ii,6);
%         clusters(ii,2) = clusters(ii,2)/clusters(ii,6);
%         clusters(ii,3) = clusters(ii,3)/clusters(ii,6);
%         clusters(ii,4) = clusters(ii,4)/clusters(ii,6);
%         clusters(ii,5) = clusters(ii,5)/clusters(ii,6);
        clusters(ii,1) = clusters(ii,1)/numpix(ii);
        clusters(ii,2) = clusters(ii,2)/numpix(ii);
        clusters(ii,3) = clusters(ii,3)/numpix(ii);
        clusters(ii,4) = clusters(ii,4)/numpix(ii);
        clusters(ii,5) = clusters(ii,5)/numpix(ii);
    end
end

% display superpixel
figure, imagesc(img);
hold on
for ii = 1:y-1
    for jj = 1:x-1
        if labels(ii,jj) ~= labels(ii,jj+1) || labels(ii,jj) ~= labels(ii+1,jj)
            plot(jj,ii,'k.');
        end
    end
end

imgVis = 0;
cIndMap = 0;

% 
time = toc;

end

