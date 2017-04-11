% function print out all front area as green 
function printFront(im, labels)

figure, imagesc(im);
hold on,
for ii = 1:size(labels,1)
    for jj = 1:size(labels,2)
        if labels(ii,jj) == 1
            plot(jj,ii,'g.');
        end
    end
end

end