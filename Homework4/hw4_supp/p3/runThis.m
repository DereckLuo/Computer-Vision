clc; close all; clear all;

% load image and display the polygon
disp('loading image and display the polygon');
im = im2double(imread('cat.jpg'));
img = imread('cat.jpg');
figure, imagesc(im);

poly = load('cat_poly.mat');

hold on,
plot(poly.poly(:,1), poly.poly(:,2), 'g*');

[y,x,z] = size(im);

%% Create current label using given poly, set labels in gch
% disp('Initialize first labeling');
% labels = zeros(y,x);
% xv = poly.poly(:,1); yv = poly.poly(:,2);
% for ii = 1:size(labels,1)
%     for jj = 1:size(labels,2)
%         if(inpolygon(jj,ii,xv,yv) == 1)
%             labels(ii,jj) = 1;
%         end
%     end
% end
% 
% save('labels.mat', 'labels');

%% Perform graphCut with given image and labels 
clc; close all;
im = im2double(imread('cat.jpg'));
labels = load('labels.mat');
labels = labels.labels;
iterations = 1;
bins = 256;

% for ii = 1:size(labels,1)
%     for jj = 1:size(labels,2)
%         if labels(ii,jj) == 1
%             im(ii,jj,1) = 0; im(ii,jj,2) = 0; im(ii,jj,3) = 225;
%         end


%     end
% end
% figure, imagesc(im);

% printFront(im, labels);

% fit GMM model for forground and background
for it = 1:iterations
    disp(it);
    [y,x] = size(labels);
    forground = zeros(y*x,3); background = zeros(y*x,3);
    fidx = 1; bidx = 1;
    for ii = 1:size(labels,1)
        for jj = 1:size(labels,2)
            if labels(ii,jj) == 1
                forground(fidx,:) = [im(ii,jj,1), im(ii,jj,2), im(ii,jj,3)];
                fidx = fidx + 1;
               
                % forground = [forground; [im(ii,jj,1), im(ii,jj,2), im(ii,jj,3), 1]];
            else 
                background(bidx, :) = [im(ii,jj,1), im(ii,jj,2), im(ii,jj,3)];
                bidx = bidx + 1;
                % background = [background; [im(ii,jj,1), im(ii,jj,2), im(ii,jj,3), 0]];
            end
        end
    end
    % remove empty rows from forground and background
    forground(all(forground==0,2),:)=[];
    background(all(background==0,2),:)=[];
    fGMM = fitgmdist(forground, 5);
    bGMM = fitgmdist(background, 5);

    % Create a new graph object 
    disp('Create new graph object');
    k1 = 1; k2 = 3; sigma = 2;
    DataCost = zeros(y,x,2);        % initialize Data cost 
    SmoothnessCost = zeros(2,2);    % initialize SoothnessCost
    vC = zeros(y,x);                % vertical neighboring cost
    hC = zeros(y,x);                % horizontal neighboring cost
    fGp = zeros(y,x);
    bGp = zeros(y,x);
    % filling DataCost matrix 

    for ii = 1:size(labels,1)
        for jj = 1:size(labels,2)
            var = [im(ii,jj,1), im(ii,jj,2), im(ii,jj,3)];
            fGp(ii,jj) = pdf(fGMM, var);
            bGp(ii,jj) = pdf(bGMM, var);
            DataCost(ii,jj,1) = -log(fGp(ii,jj)/bGp(ii,jj));
            if ii ~= y
                cur_color = var;
                next_color = [im(ii+1,jj,1),im(ii+1,jj,2),im(ii+1,jj,3)];
                vC(ii,jj) = k1 + k2*exp((-(norm(cur_color-next_color))^2)/(2*sigma^2));
            end
            if jj ~= x
                cur_color = var;
                next_color = [im(ii,jj+1,1), im(ii,jj+1,2), im(ii,jj+1,3)];
                hC(ii,jj) = k1 + k2*exp((-(norm(cur_color-next_color))^2)/(2*sigma^2));
            end
        end 
    end 

    % figure, imagesc(DataCost(:,:,1));
    % figure, imagesc(DataCost(:,:,2));

    % save('DataCost2.mat','DataCost');
    SmoothnessCost(1,2) = k1; SmoothnessCost(2,1) = k1;
    [gch] = GraphCut('open', DataCost, SmoothnessCost, vC, hC);

    % Set labels in gch 
    disp('Set labels in gch');
    [gch] = GraphCut('set', gch, labels);

    % Get current values of energy term
    [gch e] = GraphCut('energy', gch);
    disp(e);
    
    % Expand the cut 
    disp('Graph Expand Cut');
    [gch labels] = GraphCut('expand', gch);

    % Alpha Beta swap 
    % disp('Graph alpha beta swap');
    % [gch labels] = GraphCut('swap', gch, 2);
end 

% replace image pixel and printout new image 
for ii = 1:size(labels,1)
    for jj = 1:size(labels,2)
        if labels(ii,jj) == 1
            im(ii,jj,1) = 0; im(ii,jj,2) = 0; im(ii,jj,3) = 225;
        end
    end
end
figure, imagesc(im);

% Close graph and release allocated resources
disp('Graph Cut operation complete');
[gch] = GraphCut('close', gch);
