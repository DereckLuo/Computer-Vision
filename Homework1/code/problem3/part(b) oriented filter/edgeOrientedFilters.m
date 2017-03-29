function bmap = edgeOrientedFilters(im)
%call helper function and perform the non-maxima suppression, and output
%the final soft edge map

[mag, theta] = orientedFilterMagnitude(im);

bmap = nonmax(mag, theta);
%bmap = edge(mag, 'canny');

end

