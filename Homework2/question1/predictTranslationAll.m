function [newXs, newYs, outarr] = predictTranslationAll(startXs, startYs, im0, im1)
% Function to predict translation point from im0 to im1
% given initial points at startXs startYs
% output after translation location newXs newYs

newXs = startXs;
newYs = startYs;
outarr = [];

% compute the directional gradient of image 
[Ix, Iy] = imgradientxy(im0);

pointSize = size(startXs);
numPoints = pointSize(2);

% compute translation for all points 
for ii = 1:numPoints
    Px = startXs(ii);
    Py = startYs(ii);
    [u,v, out] = predictTranslation(Px, Py, Ix, Iy, im0, im1);
    if (out == true)
        outarr = [outarr, ii];
    end
    newXs(ii) =  u;
    newYs(ii) =  v;
end


end

