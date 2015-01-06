function [sal] = MotionSaliency(flow,spsMap)

magBinSz = 20;
angBinSz = 36;

imSz = size(flow(:,:,1));

disp('global bin');
tic
f_u = flow(:,:,1);
f_v = flow(:,:,2);

[f_ang, f_mag] = cart2pol(f_u, f_v);

[binFlow, e, m, binLoc] = histcn([f_mag(:) f_ang(:)], magBinSz, angBinSz);

% [bin_mag, edges_mag, mid_mag, loc_mag] = histcn(f_mag(:), 3);
% [bin_ang, edges_ang, mid_ang, loc_ang] = histcn(f_mag(:), 3);

mean_mag = histcn([f_mag(:) f_ang(:)], magBinSz, angBinSz, 'AccumData', f_mag(:), 'Fun', @mean);
mean_ang = histcn([f_mag(:) f_ang(:)], magBinSz, angBinSz, 'AccumData', f_ang(:), 'Fun', @mean);

% mean_mag = histcn(f_mag(:), 10, 'AccumData', f_mag(:), 'Fun', @mean); 
% mean_ang = histcn(f_ang(:), 16, 'AccumData', f_ang(:), 'Fun', @mean); 

% mag_table = reshape(mean_mag(loc_mag), size(f_mag));
% ang_table = reshape(mean_ang(loc_ang), size(f_ang));

binnedMeanMag = mean_mag(sub2ind([magBinSz,angBinSz], binLoc(:,1), binLoc(:,2)));
figure;imagesc(reshape(binnedMeanMag, size(f_u)));
% binnedMeanAng = mean_ang(sub2ind([magBinSz,angBinSz], binLoc(:,1), binLoc(:,2)));

[uMeanFlow, vMeanFlow] = pol2cart(mean_ang(:), mean_mag(:)); % careful, theta first arg
meanFlow = [uMeanFlow, vMeanFlow];
          
normBinFlow = binFlow/sum(binFlow(:)); % normalize to 1
toc

disp('sps gen');
tic;
spsInds = unique(spsMap);
spsInds = spsInds(spsInds > 0);
totSpsCount = length(spsInds);
sps(totSpsCount) = struct('pixelCount',[],'meanFlow',[],'center',[],'binCount',[],'index',[]);
spsCount = 1;

for i = spsInds'
    spsLoc = spsMap == i;
    
    inds = find(spsLoc);
    [y,x] = ind2sub(imSz, inds);
    sps(spsCount).center = [sum(y)/length(inds) sum(x)/length(inds)];
    
    spsBinnedFlowVals = binLoc(inds,:);
    
    spsFlowVals = [mean_mag(sub2ind([magBinSz,angBinSz], spsBinnedFlowVals(:,1), spsBinnedFlowVals(:,2))) ...
                  mean_ang(sub2ind([magBinSz,angBinSz], spsBinnedFlowVals(:,1), spsBinnedFlowVals(:,2)))];
              
    [sps_u_norm, sps_v_norm] = pol2cart(spsFlowVals(:,2), spsFlowVals(:,1));
              
    spsBinCount = zeros(magBinSz,angBinSz);
    spsBinInds = sub2ind([magBinSz,angBinSz], spsBinnedFlowVals(:,1), spsBinnedFlowVals(:,2));
    uSpsBinInds = unique(spsBinInds);
    counts = arrayfun(@(y) length(spsBinInds(spsBinInds == y)), uSpsBinInds);
    spsBinCount(uSpsBinInds) = counts;
    
    sps(spsCount).pixelCount = nnz(spsLoc);
%     sps(spsCount).meanFlow = mean(spsFlowVals,1);
    sps(spsCount).meanFlow = mean([sps_u_norm, sps_v_norm],1);
    sps(spsCount).binCount = spsBinCount/sum(spsBinCount(:));
    sps(spsCount).index = i;
    spsCount = spsCount + 1;
end
toc;

disp('score gen');
tic
S_md = zeros(totSpsCount,1);
for i = 1:totSpsCount
%     tic
    u1 = repmat(meanFlow(:,1), 1, magBinSz*angBinSz);
    u2 = repmat(meanFlow(:,1)', magBinSz*angBinSz, 1);
    v1 = repmat(meanFlow(:,2), 1, magBinSz*angBinSz);
    v2 = repmat(meanFlow(:,2)', magBinSz*angBinSz, 1);
    u3 = u1 - u2;
    v3 = v1 - v2;
    uvNorm = sqrt(u3.^2 + v3.^2);
    nBF1 = repmat(normBinFlow(:), 1, magBinSz*angBinSz);
    inSum = sum(uvNorm.*nBF1,1);
    outSum = sum(sps(i).binCount(:)'.*inSum);
%     toc
    
    tic
    outerSum = 0;
    for j = 1:magBinSz*angBinSz
        
        innerSum = 0;
        for k = 1:magBinSz*angBinSz
            innerSum = innerSum + norm(meanFlow(j,:) - meanFlow(k,:)) * normBinFlow(k);
        end
        outerSum = outerSum + sps(i).binCount(j) * innerSum;
    end
    toc
    
    S_md(i) = outSum;
end
toc

% tmp = zeros(size(spsMap));
% for i=1:totSpsCount
%     spsLoc = spsMap == sps(i).index;
%     tmp(spsLoc) = S_md(i);
% end
% figure;imagesc(tmp);

disp('adjust score');


end 