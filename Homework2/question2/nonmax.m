function [keyXs, keyYs] = nonmax(im, tau)
% nonmax suppression over image, chose point over threshold tau
% and largest in sub 5x5 block.

keyXs = [];
keyYs = [];
[row, column] = size(im);

for ii = 3:row-2
    for jj = 3:column-2
        current = im(ii,jj);
        if(current >= tau)
            local = [im(ii-2,jj-2),im(ii-2,jj-1),im(ii-2,jj),im(ii-2,jj+1),im(ii-2,jj+2),im(ii-1,jj-2),im(ii-1,jj-1),im(ii-1,jj),im(ii-1,jj+1),im(ii-1,jj+2),im(ii  ,jj-2),im(ii  ,jj-1),im(ii  ,jj+1),im(ii  ,jj+2),im(ii+1,jj-2),im(ii+1,jj-1),im(ii+1,jj),im(ii+1,jj+1),im(ii+1,jj+2),im(ii+2,jj-2),im(ii+2,jj-1),im(ii+2,jj),im(ii+2,jj+1),im(ii+2,jj+2)];
            local_max = max(local);
            if(current >= local_max)
                keyXs = [keyXs, jj];
                keyYs = [keyYs, ii];
            end
        end 
    end 
end 


end

