function neighbouringSpInds = findNeighbouringSuperpixelIndices(spIndex, frame, spArray)
%findNeighbouringSuperpixelIndices Takes a superpixel index and frame containing an array of
%superpixels and returns the neighbouring superpixel indices

if isempty(spArray)
    spArray = frame.spArray;
end

% % use the x and y inds of the closest superpixel to create a bounding box around it
% minX = max(min(spArray(spIndex).xInds) - 1, 1);
% maxX = min(max(spArray(spIndex).xInds) + 1, frame.imgSize(2));
% minY = max(min(spArray(spIndex).yInds) - 1, 1);
% maxY = min(max(spArray(spIndex).yInds) + 1, frame.imgSize(1));
% 
% % get the indices of the neighbouring superpixels in spArray
% neighbouringSpInds = find(ismember(frame.spInds, unique(frame.spMap(minY:maxY, minX:maxX))) > 0);

[p(:,:,2), p(:,:,1)] = meshgrid(-2:2,-2:2);
p = reshape(p,[],2);
p = p([1:12 14:25], :);

glcms = graycomatrix(frame.spMap, 'NumLevels', frame.spNum+1, 'GrayLimits', [0 frame.spNum], 'Offset', p);

end

