% MATLAB implementation of the color frequency image.
% Code inspired by:
%-------------------------------------------------------------------------------
% http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?tp=&arnumber=4421050&isnumber=4420928
% Introduction of frequency image and applications
% Kashiwagi, T.   Oe, S.   
% Tokushima Prefectural Ind. Technol. Center, Tokushima;
% 
% This paper appears in: SICE (Society of Instrument and Control Engineers), 2007 Annual Conference
% Publication Date: 17-20 Sept. 2007
% On page(s): 584-591
% Location: Takamatsu, 
% ISBN: 978-4-907764-27-2
% INSPEC Accession Number: 9853526
% Digital Object Identifier: 10.1109/SICE.2007.4421050
% Current Version Published: 2008-01-07 
%  
% Abstract
% We introduce a new imaging method of "Frequency Image" and its applications. 
% This new image is made from multi-dimensional features of an image and presents 
% a distribution of the frequency of feature vectors. Especially, selecting R, G and B colors 
% as the three features of an image, Frequency Image is a feature image of which pixel 
% represents the frequency of the same color pixels. In this paper, first, we explain 
% the basic idea and how to make this Frequency Image. Next, we introduce some effective 
% applications using this image. Finally we mention excellent potentials of Frequency Images. 
%-------------------------------------------------------------------------------
%
% This MATLAB implementation is written by Mark Hayworth, Ph.D.
% Only tested on 24 bit (8 bit RGB) images.
% Imaging Core Section
% email: hayworth dot ms at pg dot com
% The Procter & Gamble Company, Cincinnati, Ohio, USA.
% June 2009.

clear all;
close all;
clc;

% Change the current folder to the folder of this m-file.
% (The line of code below is from Brett Shoelson of The Mathworks.)
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end
	
% For convenience in browsing, set a starting folder from which to browse.
startingFolder = 'C:\Documents and Settings\tk2013\My Documents\My Pictures\Macbeth';
if ~exist(startingFolder, 'dir')
	% If that folder doesn't exist, just start in the current folder.
	startingFolder = pwd;
end
% Get the name of the color image file that the user wants to use.
defaultFileName = fullfile(startingFolder, '*.*');
[baseFileName, folder, FilterIndex] = uigetfile(defaultFileName, 'Select a color image file');
if baseFileName == 0
	return;
end

imageFullFileName = fullfile(folder, baseFileName);

% imageFullFileName = 'C:\Documents and Settings\tk2013\My Documents\My Pictures\Macbeth\Macbeth blurry noisy.tif';
originalImage = imread(imageFullFileName);
% Determine it's size and turn it into color if it's not already color.
[rows columns numberOfColorBands] = size(originalImage);
if numberOfColorBands < 3
	% It's not a color image - it's a monochrome image.  Turn it into a color image.
	uiwait(msgbox('Note: this is a monochrome image, not a color image.'));
	originalImage = cat(3, originalImage, originalImage, originalImage);
end
% Display the original color image.
hFigOriginal = figure;
set(hFigOriginal, 'Name', 'originalImage');
imshow(originalImage, []);
fontSize = 18;
title('Original Color Image', 'FontSize', fontSize);
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%drawnow;

% Construct the 3D gamut.
lutSize = 256;  % Use 256 to get maximum resolution possible out of a 24 bit RGB image.
lutSize = 8;  % Use a smaller LUT size to get colors grouped into fewer, larger groups in the gamut.
reductionFactor = double(256) / double(lutSize);
gamut3D = zeros(lutSize, lutSize, lutSize);
for row = 1 : rows
	for col = 1: columns
		redValue = floor(double(originalImage(row, col, 1)) / reductionFactor) + 1;	% Convert from 0-255 to 1-256
		greenValue = floor(double(originalImage(row, col, 2)) / reductionFactor) + 1;	% Convert from 0-255 to 1-256
		blueValue = floor(double(originalImage(row, col, 3)) / reductionFactor) + 1;	% Convert from 0-255 to 1-256
		gamut3D(redValue, greenValue, blueValue) = gamut3D(redValue, greenValue, blueValue) + 1;
	end
end

% Now construct the color frequency image.
% Make an image where we get the color of the original image, and have the output value of the color
% frequency image be the number of times that exact color occurred in the original image.
colorFrequencyImage = zeros(rows, columns);
for row = 1 : rows
	for col = 1: columns
		redValue = floor(double(originalImage(row, col, 1)) / reductionFactor) + 1;	% Convert from 0-255 to 1-256
		greenValue = floor(double(originalImage(row, col, 2)) / reductionFactor) + 1;	% Convert from 0-255 to 1-256
		blueValue = floor(double(originalImage(row, col, 3)) / reductionFactor) + 1;	% Convert from 0-255 to 1-256
		freq = gamut3D(redValue, greenValue, blueValue);
		colorFrequencyImage(row, col) =  freq;
	end
end
numberOfColors = double(rows) * double(columns);
disp(numberOfColors);
numberOfColors = sum(sum(sum(gamut3D)));
disp(numberOfColors);
sumr = sum(sum(colorFrequencyImage));
disp(sumr);
figure;
imshow(colorFrequencyImage, []);
normImage = colorFrequencyImage * 255 / max(max(colorFrequencyImage));
% imtool(colorFrequencyImage, [0 25]);
title('Color Frequency Image', 'FontSize', fontSize);
%imwrite(normImage, 'Color Frequency Output Image.tif');
% [lowThreshold highThreshold lastThresholdedBand] = threshold(1, 100, colorFrequencyImage);
message = sprintf('Done!\nThe color frequency image represents how often\nthe color of each pixel occurs anywhere in the image.\nIf the pixel is brighter then that color appears more often in the image.');
uiwait(msgbox(message));

