function [spArray, motionSignificanceHistory, saliencySpatial] = populateSpArray(frameNow, framePrev)
%POPULATESPARRAY Populates the superpixel array by looping through each superpixel and computing its
%histograms and properties

% motion significance window length
motionSignigicanceWinLength = 5;
motionSignificanceHistory = [];
pixelSaliencySpatial = zeros(frameNow.imgSize);

% preallocate spArray
spArray = repmat(struct( ...
    'pixelCount',[], ...
    'yInds', [], ...
    'xInds', [], ...
    'neighbourSpInds', [], ...
    'center',[], ...
    'pixelVals',[], ...
    'flowVals',[], ...
    'motionHist',[], ...
    'colorHist',[], ...
    'meanMotion',[], ...
    'interFrameSimilarity', [], ...
    'saliencyMotionDistinctiveness', [], ...
    'saliencyTemporalPrediction', [], ...
    'saliencyTemporal', [], ...
    'prevMatchingSpInd', [], ...
    'intraFrameSimilarity', zeros(frameNow.spNum, 1), ...
    'saliencyGlobalContrast', [], ...
    'saliencySpatialSparcity', [], ...
    'saliencySpatial', []), frameNow.spNum, 1);

if frameNow.frameNum ~= 1
    prevCenters = reshape([framePrev.spArray.center], 2, [])';
end

[offsets(:,:,2), offsets(:,:,1)] = meshgrid(-2:2,-2:2);
offsets = reshape(offsets,[],2);
offsets = offsets([1:6 10 11 15 16 20:25], :); % remove unnecessary center offsets (5x5 window, only need outer edge due to 0 border pixels)
spGrayCoMatrix = graycomatrix(frameNow.spMap, 'NumLevels', frameNow.spNum+1, 'GrayLimits', [0 frameNow.spNum], 'Offset', offsets);
spGrayCoMatrix = sum(spGrayCoMatrix(2:end, 2:end, :), 3);

% loop through each superpixel
for i = 1:frameNow.spNum
    % get indices and count
    spArray(i).pixelInds = find(frameNow.spMap == frameNow.spInds(i));
    spArray(i).pixelCount = length(spArray(i).pixelInds);
    
    % compute neighbours
    spArray(i).neighbourSpInds = find(spGrayCoMatrix(:,i) > 0);
    %spArray(i).neighbourSpInds = spArray(i).neighbourSpInds(spArray(i).neighbourSpInds ~= i);
    
    % compute center
    [spArray(i).yInds, spArray(i).xInds] = ind2sub(frameNow.imgSize, spArray(i).pixelInds);
    spArray(i).center = [sum(spArray(i).yInds)/length(spArray(i).pixelInds) sum(spArray(i).xInds)/length(spArray(i).pixelInds)];
    spArray(i).distanceToFrameCenter = pdist2(spArray(i).center, frameNow.center, 'euclidean');
    
    % compute superpixel level color histogram and normalize
    [spArray(i).colorHist, ~, ~, spArray(i).colorLoc] = histcn(frameNow.labImg(spArray(i).pixelInds,:), frameNow.colorEdges{1}, frameNow.colorEdges{2}, frameNow.colorEdges{3});
    spArray(i).colorHist = spArray(i).colorHist/sum(spArray(i).colorHist(:));
    
    % compute superpixel level motion histogram and normalize
    if frameNow.frameNum ~= 1
        spArray(i).motionHist = histcn(frameNow.flowPolar(spArray(i).pixelInds,:), frameNow.motionEdgesPolar{1}, frameNow.motionEdgesPolar{2});
        spArray(i).motionHist = spArray(i).motionHist/sum(spArray(i).motionHist(:));
        
        % calculate mean flow of the superpixel
        spMeanFlow = mean(frameNow.flowCart(spArray(i).pixelInds,:), 1);
        % project the center backwards
        projectedCenter = repmat(spArray(i).center - spMeanFlow, framePrev.spNum, 1);
        % compute distances between projected center and all previous superpixel centers
        norms = sqrt(sum((projectedCenter - prevCenters).^2, 2));
        % find the closest match
        [~, spArray(i).prevMatchingSpInd] = min(norms);
        % find neighbouring indices of the most likely superpixel from the previous frame
        %spArray(i).prevNeighbouringSpInds = findNeighbouringSuperpixelIndices(minSpIndex, framePrev, []);
        % calculate interframe similarity between current superpixel and each neighbour
        
        prevNeighbouringSpInds = framePrev.spArray(spArray(i).prevMatchingSpInd).neighbourSpInds;
        spArray(i).interFrameSimilarity = zeros(length(prevNeighbouringSpInds), 1);
        currSp = spArray(i);
        for j = 1:length(prevNeighbouringSpInds)
            prevSp = framePrev.spArray(prevNeighbouringSpInds(j));
            spArray(i).interFrameSimilarity(j) = sum(sqrt(currSp.colorHist(:) .* prevSp.colorHist(:))) * ...
                (1 - norm(currSp.center - spMeanFlow - prevSp.center)/framePrev.diagonal);
        end
    end
end

distancesToFrameCenter = [spArray.distanceToFrameCenter]';

% frameNow.spNum ~ 700
% frameNow.contrastBinCount = 16 (per channel)
% frameNow.colorMean = 4096x3 (mean color for each bin)
% frameNow.colorHist = 16x16x16 (global histogram counts)
% spArray(i).colorHist = 16x16x16 (superpixel histogram counts)

% arrange the colorHists from 16x16x16 into a mat of 4096xspNum
colorHistsMat = cell2mat(cellfun(@(x) x(:), {spArray.colorHist}, 'UniformOutput', false));

% EQ 8
% compute the static global weighted distances
colorDistanceMatrix = pdist2(frameNow.colorMean, frameNow.colorMean, 'euclidean');
weightedColorDistanceMatrix = colorDistanceMatrix * frameNow.colorHist(:);

% compute the global contrast saliency and assign back into spArray
saliencyGlobalContrast = colorHistsMat'*weightedColorDistanceMatrix;

% EQ 9
% compute the bhattacharyya coefficient
colorHistsSqrt = sqrt(colorHistsMat);
bhattacharyyaCoeff = colorHistsSqrt'*colorHistsSqrt;

% compute the weights based on distances from centers
centers = reshape([spArray.center], 2, [])';
weightedCentersDistanceMatrix = 1 - pdist2(centers, centers, 'euclidean')/frameNow.diagonal;

% EQ 10
% compute the intra frame similarity and the color spatial distribution
intraFrameSimilarity = bhattacharyyaCoeff.*weightedCentersDistanceMatrix;
colorSpatialDistribution = (intraFrameSimilarity*distancesToFrameCenter)./sum(intraFrameSimilarity, 2);

% EQ 11 & 12
% calculate spatial sparcity saliency and spatial saliency
maxColorSpatialDistribution = max(colorSpatialDistribution);
minColorSpatialDistribution = min(colorSpatialDistribution);

saliencySpatialSparcity = (maxColorSpatialDistribution - colorSpatialDistribution)/(maxColorSpatialDistribution - minColorSpatialDistribution);
saliencySpatial = num2cell(saliencySpatialSparcity.*saliencyGlobalContrast);
[spArray(:).saliencySpatial] = saliencySpatial{:};

% TODO: check this, feels wrong
for i = 1:frameNow.spNum
    % find neighbour indices for each superpixel in the current frame
    %spArray(i).nowNeighbouringSpInds = findNeighbouringSuperpixelIndices(i, frameNow, spArray);
    % get the bin indices of each pixel in the superpixel
    colorBinInds = sub2ind(repmat(frameNow.contrastBinCount, 1, 3), spArray(i).colorLoc(:,1), spArray(i).colorLoc(:,2), spArray(i).colorLoc(:,3));
    % compute the pixel-level spatial saliency score for each pixel in the superpixel
    numeratorSum = zeros(spArray(i).pixelCount, 1);
    denomenatorSum = zeros(spArray(i).pixelCount, 1);
    for j = 1:length(spArray(i).neighbourSpInds)
        numeratorSum = numeratorSum + (spArray(j).saliencySpatial * spArray(j).colorHist(colorBinInds));
        denomenatorSum = denomenatorSum + spArray(j).colorHist(colorBinInds);
    end
    pixelSaliencySpatial(spArray(i).pixelInds) = numeratorSum./denomenatorSum;
end


if frameNow.frameNum ~= 1
    % calculate motion distinctiveness
    % vectorized form of the motion distinctiveness calculation. It flattens out the double loop
    % into large q^2 x q^2 matrices, computes the inner sum for all combinations of j and k,
    % multiplies each value by the correspoing global motion bin count at k, sums that value, and
    % then does the outer j sum and multiply with the superpixel motion bin count. Higher memory
    % requirement but much, much faster.
    
    motionHistsMat = cell2mat(cellfun(@(x) x(:), {spArray.motionHist}, 'UniformOutput', false));
    
    flowDistanceMatrix = pdist2(frameNow.motionMeanCart, frameNow.motionMeanCart, 'euclidean');
    weightedFlowDistanceMatrix = flowDistanceMatrix * frameNow.motionHist(:);
    
    saliencyMotionDistinctiveness = motionHistsMat'*weightedFlowDistanceMatrix;
    
    
    % calculate motion significance
    motionSignificance = [spArray.pixelCount] * saliencyMotionDistinctiveness;
    if frameNow.frameNum == 2
        motionSignificanceHistory = motionSignificance;
    else
        % calculate the adjustment weight based on motion significance for each frame
        adjustmentWeight = max(median(framePrev.motionSignificanceHistory)/(median(framePrev.motionSignificanceHistory) + motionSignificance), 0.5);
        motionSignificanceHistory = [framePrev.motionSignificanceHistory motionSignificance];
        motionSignificanceHistory = motionSignificanceHistory(max(length(motionSignificanceHistory)- motionSignigicanceWinLength, 1):end);
    end
    
    % TODO: double loop, see if I can't collapse it/vectorize it
    
    
    
    if frameNow.frameNum == 2
        salMotionDistCell = num2cell(saliencyMotionDistinctiveness);
        [spArray(:).saliencyTemporal] = salMotionDistCell{:};
    else
        prevMatchingSpInds = [spArray.prevMatchingSpInd];
        prevNeighbouringSpInds = {framePrev.spArray(prevMatchingSpInds).neighbourSpInds};
        interFrameSimilarity = {spArray.interFrameSimilarity};

        saliencyTemporalPrediction = cellfun(@(inds, interFrameSim) ([framePrev.spArray(inds).saliencyTemporal]*interFrameSim)/sum(interFrameSim), prevNeighbouringSpInds, interFrameSimilarity, 'UniformOutput', false);
        emptyInds = cellfun('isempty', saliencyTemporalPrediction);
        nanInds = cellfun(@(x) any(isnan(x)), saliencyTemporalPrediction);
        saliencyTemporalPrediction(emptyInds) = {0};
        saliencyTemporalPrediction(nanInds) = {0};
        saliencyTemporalPrediction = cell2mat(saliencyTemporalPrediction);
        
        saliencyTemporal = num2cell(adjustmentWeight*saliencyTemporalPrediction' + ...
                    (1 - adjustmentWeight)*saliencyMotionDistinctiveness);
        [spArray(:).saliencyTemporal] = saliencyTemporal{:};
    end
end

end
