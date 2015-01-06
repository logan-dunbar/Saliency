function [binColors, binLoc, meanColorsL, meanColorsA, meanColorsB, binnedImg] = GlobalColorHist(img, q, alpha)
%GLOBALCOLORHIST Take an CIEL*A*B* image and compute the global color
%histogram
%   Takes a CIEL*a*b* image and quantizes the color histogram into qxqxq bins, while
%   keeping the top alpha percent of high probability colors.


l_bin_edges = 0:100/q:100;
ab_bin_edges = -128:255/q:127;

imgReshape = reshape(img, [], 3);

% [count, edges, ~, loc] = histcn(imgReshape, l_bin_edges, ab_bin_edges, ab_bin_edges);

% tic
% binPerm = PermsRep(1:q,3);
% binColors = struct([]);
% for i=1:q^3
%     locInds = ismember(loc, binPerm(i,:),'rows');
%     binColors(i).meanColor = mean(imgReshape(locInds,:),1);
% end
% toc

% tic
% binColors = struct([]);
% for i=1:q
%     for j=1:q
%         for k=1:q
%             locInds = [loc(:,1) == i loc(:,2) == j loc(:,3) == k];
%             locInds = all(locInds>0, 2);
%             vals = imgReshape(locInds,:);
%             binColors(i,j,k).meanColor = mean(vals,1);
%             binColors(i,j,k).index = [i j k];
%             binColors(i,j,k).edges = [edges{1}(i:i+1)' edges{2}(j:j+1)' edges{3}(k:k+1)'];
%             binColors(i,j,k).count = count(i,j,k);
%         end
%     end
% end
% toc

tic
[binColors, ~, ~, binLoc] = histcn(imgReshape, l_bin_edges, ab_bin_edges, ab_bin_edges);
meanColorsL = histcn(imgReshape, l_bin_edges, ab_bin_edges, ab_bin_edges, 'AccumData', imgReshape(:,1), 'Fun', @mean);
meanColorsA = histcn(imgReshape, l_bin_edges, ab_bin_edges, ab_bin_edges, 'AccumData', imgReshape(:,2), 'Fun', @mean);
meanColorsB = histcn(imgReshape, l_bin_edges, ab_bin_edges, ab_bin_edges, 'AccumData', imgReshape(:,3), 'Fun', @mean);

binnedImg = [meanColorsL(sub2ind([q,q,q], binLoc(:,1), binLoc(:,2), binLoc(:,3))) ...
             meanColorsA(sub2ind([q,q,q], binLoc(:,1), binLoc(:,2), binLoc(:,3))) ...
             meanColorsB(sub2ind([q,q,q], binLoc(:,1), binLoc(:,2), binLoc(:,3)))];
toc
             
% l = img(:,:,1);
% a = img(:,:,2);
% b = img(:,:,3);
% 
% [~,l_bin] = histc(l(:),l_bin_edges);
% [~,a_bin] = histc(a(:),ab_bin_edges);
% [~,b_bin] = histc(b(:),ab_bin_edges);
% 
% glMap.l_bins = struct([]);
% glMap.a_bins = struct([]);
% glMap.b_bins = struct([]);
% for i = 1:q
%     glMap.l_bins(i).vals = l(l_bin == i);
%     glMap.l_bins(i).count = length(glMap.l_bins(i).vals);
%     glMap.l_bins(i).mean_color = mean(glMap.l_bins(i).vals);
%     glMap.a_bins(i).vals = l(l_bin == i);
%     glMap.a_bins(i).count = length(glMap.a_bins(i).vals);
%     glMap.a_bins(i).mean_color = mean(glMap.a_bins(i).vals);
%     glMap.b_bins(i).vals = l(l_bin == i);
%     glMap.b_bins(i).count = length(glMap.b_bins(i).vals);
%     glMap.b_bins(i).mean_color = mean(glMap.b_bins(i).vals);
% end

% [inds,inds] = sort(
% glHist

end