function [colorHist, colorEdges, colorMean] = populateContrast(frame)
%POPULATEGLOBALCONTRASTHIST Populates the contrast histogram, quantized contrast table and pixel bin
%locations for an image

labImg = frame.labImg;
binCount = frame.contrastBinCount;

% compute histogram and normalize
[colorHist, colorEdges, ~, ~] = histcn(labImg, binCount, binCount, binCount);
colorHist = colorHist/sum(colorHist(:));

% compute mean color for each bin
meanL = histcn(labImg, binCount, binCount, binCount, 'AccumData', labImg(:,1), 'Fun', @mean);
meanA = histcn(labImg, binCount, binCount, binCount, 'AccumData', labImg(:,2), 'Fun', @mean);
meanB = histcn(labImg, binCount, binCount, binCount, 'AccumData', labImg(:,3), 'Fun', @mean);
colorMean = [meanL(:) meanA(:) meanB(:)];

end

