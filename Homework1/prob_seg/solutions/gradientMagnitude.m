function [mag, theta] = gradientMagnitude(im, sigma)

% step 1 : smooth the image with Gaussian std = sigma
smooth = imgaussfilt(im, sigma);

[Rmag, Rdir] = imgradient(smooth(:,:,1));
[Gmag, Gdir] = imgradient(smooth(:,:,2));
[Bmag, Bdir] = imgradient(smooth(:,:,3));
[Imag, Idir] = imgradient(rgb2gray(smooth));

imSize = size(Gmag);
mag = zeros(imSize);
theta = zeros(imSize);

row = imSize(1);
colum = imSize(2);

for ii = 1:row
    for jj = 1:colum
        cell = [Rmag(ii,jj), Gmag(ii,jj), Bmag(ii,jj), Imag(ii,jj)];
        maxmag = max(cell);
        if maxmag == cell(1)
            theta(ii,jj) = Rdir(ii,jj);
        elseif maxmag == cell(2)
            theta(ii,jj) = Gdir(ii,jj);
        elseif maxmag == cell(3)
            theta(ii,jj) = Bdir(ii,jj);
        else
            theta(ii,jj) = Idir(ii,jj);
        end
        mag(ii,jj) = norm(cell);
    end
end
