function [newPx, newPy, out] = predictTranslation(Px, Py, Ix, Iy, im0, im1)
% predict translation on a single point input

newPx = Px; newPy = Py; out = false;
% calculate Ix, and Iy using interperlation 
window_size = -7:7;
[X,Y] = meshgrid(window_size, window_size);
X0 = X + Px; Y0 = Y + Py;
WIx = interp2(Ix, X0, Y0);
WIy = interp2(Iy, X0, Y0);
if(isnan(WIx(1,1)) || isnan(WIx(1,15)) || isnan(WIx(15,1)) || isnan(WIx(15,15)))
    out = true;
    return;
end
% forming matrix A
WIx = reshape(WIx, 15*15, 1);
WIy = reshape(WIy, 15*15, 1);
WIx = WIx'; WIx = WIx(:);
WIy = WIy'; WIy = WIy(:);
A = [WIx, WIy];

% calculate initial It0 and It1
It0 = interp2(im0, X0, Y0);
It1 = interp2(im1, X0, Y0);

% define x and y displacement
u = 0; v = 0;

% perform 20 iterations of prediction transfered point
for ii = 1:25
    % compute It 
    b = It1 - It0;
    if(isnan(b(1,1)) || isnan(b(1,15)) || isnan(b(15,1)) || isnan(b(15,15)))
        out = true;
        return;
    end 
    b = reshape(b, 15*15, 1);
    b = b'; b = b(:);
    AtA = (A')*A; 
    Atb = (A')*b;
    
    displacement = linsolve(AtA, (-1*Atb));
    %displacement = AtA\(-1*Atb);
    u = displacement(1); v = displacement(2);
    % compute predicted location 
    newPx = newPx + u; newPy = newPy + v; % each prediction is a correctness added to the origional prediction 
    % compute new It1
    X0 = X + newPx; Y0 = Y + newPy;
    It1 = interp2(im1, X0, Y0);
end

