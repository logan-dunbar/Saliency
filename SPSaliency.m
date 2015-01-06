function salImg = SPSaliency(img, spsMap,num)
%SPSALIENCY Compute a saliency map based on superpixels
%   Takes as input an image and its custered superpixels and computes a
%   saliency map based on color and distance similarities

% disp('init');
% tic;
q = 16;

img = double(img);
imSz = size(img);
numSps = max(spsMap(:));
spsCount = 1;

labImg = RGB2Lab(img);
% binnedLabImg(:,:,1) = ceil(labImg(:,:,1)/(100/q));
% binnedLabImg(:,:,2) = ceil((labImg(:,:,2) + 128)/(255/q));
% binnedLabImg(:,:,3) = ceil((labImg(:,:,3) + 128)/(255/q));
% binnedLabImg(binnedLabImg == 0) = 1; % fix any 0 entries (ceil(0) = 0)
% binnedLabImg = reshape(binnedLabImg,[],3);
% toc;

disp('hist gen');
tic;
[binColors, binLoc, meanColorsL, meanColorsA, meanColorsB, binnedImg] = GlobalColorHist(labImg, q, 0.95);
% activeBins = binColors([binColors.count] > 0);
% [~, sortedInds] = sort([activeBins.count], 'descend');
% sortedBinColors = activeBins(sortedInds);
toc

%     imRs = reshape(binnedLabImg,[],3);
%     for i = 1:imSz(1)*imSz(2)
%         myImg(i,:) = binColors(imRs(i,1), imRs(i,2), imRs(i,3)).meanColor;
%     end
%     myImgRs = reshape(myImg,imSz);
%     figure;imshow(Lab2RGB(myImgRs));

disp('sps gen');
tic;
spsInds = unique(spsMap);
spsInds = spsInds(spsInds > 0);
totSpsCount = length(spsInds);
sps(totSpsCount) = struct('pixelCount',[],'meanColor',[],'center',[],'binCount',[],'index',[]);
spsCount = 1;

% mySpsPixCount = zeros(totSpsCount,1);
% mySpsMeanColor = zeros(totSpsCount,3);
% mySpsBinCount = zeros(totSpsCount,q*q*q);
% mySpsIndex = zeros(totSpsCount,1);

for i = spsInds'
    spsLoc = spsMap == i;
    
    %sneaky
    %builtin('_brace', num2cell(magic(5)), 3, 3)
    %builtin('_paren', img(:,:,1), spsBin);
    
    inds = find(spsLoc);
    [y,x] = ind2sub(imSz, inds);
    sps(spsCount).center = [sum(y)/length(inds) sum(x)/length(inds)];
    
    spsBinnedPixVals = binLoc(inds,:);
    
    spsPixVals = [meanColorsL(sub2ind([q,q,q], spsBinnedPixVals(:,1), spsBinnedPixVals(:,2), spsBinnedPixVals(:,3))) ...
                   meanColorsA(sub2ind([q,q,q], spsBinnedPixVals(:,1), spsBinnedPixVals(:,2), spsBinnedPixVals(:,3))) ...
                   meanColorsB(sub2ind([q,q,q], spsBinnedPixVals(:,1), spsBinnedPixVals(:,2), spsBinnedPixVals(:,3)))];
    
    spsBinCount = zeros(q,q,q);
    spsBinInds = sub2ind([q,q,q], spsBinnedPixVals(:,1), spsBinnedPixVals(:,2), spsBinnedPixVals(:,3));
    uSpsBinInds = unique(spsBinInds);
    counts = arrayfun(@(y) length(spsBinInds(spsBinInds == y)), uSpsBinInds);
    spsBinCount(uSpsBinInds) = counts;
    
    
%     spsPixVals = zeros(size(spsBinnedPixVals));
%     spsBinCount = zeros(q,q,q);
%     for j = 1:size(spsBinnedPixVals,1)
%         spsPixVals(j,:) = binColors(spsBinnedPixVals(j,1), spsBinnedPixVals(j,2), spsBinnedPixVals(j,3)).meanColor;
%         spsBinCount(spsBinnedPixVals(j,1), spsBinnedPixVals(j,2), spsBinnedPixVals(j,3)) = spsBinCount(spsBinnedPixVals(j,1), spsBinnedPixVals(j,2), spsBinnedPixVals(j,3)) + 1;
%     end
    

%     mySpsPixCount(spsCount) = nnz(spsLoc);
%     mySpsMeanColor(spsCount,:) = mean(spsPixVals,1);
%     normSpsBinCount = spsBinCount/sum(spsBinCount(:));
%     mySpsBinCount(spsCount,:) = normSpsBinCount(:);
%     mySpsIndex(spsCount) = i;

    sps(spsCount).pixelCount = nnz(spsLoc);
    sps(spsCount).meanColor = mean(spsPixVals,1);
    sps(spsCount).binCount = spsBinCount/sum(spsBinCount(:));
    sps(spsCount).index = i;
    spsCount = spsCount + 1;
end
toc



disp('score gen');
tic;
simc = zeros(totSpsCount, totSpsCount);
simd = zeros(totSpsCount, totSpsCount);
W = zeros(totSpsCount, totSpsCount);
GC = zeros(totSpsCount,1);
SS = zeros(totSpsCount,1);
d = norm(imSz(1:2));
imCenter = imSz(1:2)/2;
for i = 1:totSpsCount
    ssNum = 0;
    ssDenom = 0;
    for j = 1:totSpsCount
        countIntersect = min(sps(i).binCount, sps(j).binCount);
        simc(i,j) = sum(countIntersect(:));
        
        simd(i,j) = 1 - norm(sps(i).center - sps(j).center)/d;
        W(i,j) = sps(j).pixelCount * simd(i,j);
        GC(i) = GC(i) + W(i,j)*norm(sps(i).meanColor - sps(j).meanColor);
        ssNum = ssNum + norm(sps(j).center - imCenter)*simc(i,j)*simd(i,j);
        ssDenom = ssDenom + simc(i,j)*simd(i,j);
    end
    SS(i) = ssNum/ssDenom;
end

sim = simc.*simd;
NGC = (GC - min(GC))/(max(GC) - min(GC));
NSS = (SS - max(SS))/(min(SS) - max(SS));

RGC = sum(sim.*repmat(NGC,1,totSpsCount), 2)./sum(sim,2);
RSS = sum(sim.*repmat(NSS,1,totSpsCount), 2)./sum(sim,2);

sal = RGC.*RSS;
toc;

disp('sal map');
tic;
salImg = zeros(imSz(1:2));
for i = 1:totSpsCount
    spsLoc = spsMap == sps(i).index;
    salImg(spsLoc) = sal(i);
end
toc;
figure;imagesc(salImg);
colormap gray;
title(num2str(num));
end