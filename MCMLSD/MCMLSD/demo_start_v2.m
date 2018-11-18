clear
clc
close all

addpath(genpath('code/'));
addpath('Imgs/');
img = imread(['Imgs', filesep, 'P1040823hr.jpg']);
img = imresize(img, [round(size(img,1)/4), round(size(img,2)/4)]);

%compute the kernel for the image size
%you only need to compute the kernal once for one an image size
[kernels, kernels_flip, kernel_params] =kernelInitialization(img);
ticId = tic;
%the lines variable contains the detected line segmentations it arranged as
%[x1 y1 x2 y2 probability]
%The fullLines are the detected lines. It is arranged as [rho theta probability]
[lines, fullLines] =lineSegmentation_HighRes(img,kernels, kernels_flip, kernel_params);
display('Total time');
toc(ticId)
fig = figure;
imshow(img);
hold all
   %Order lines by probability
   lines = sortrows(lines, -5);
   ttlLines = size(lines,1);
   for i = 1:90
     %plot the top 90 lines
     line([lines(i,1) lines(i,3)], [lines(i,2) lines(i,4)],'Color', rand(1,3), 'LineWidth', 3);
   end

%please use code in Evaluation code.zip to evaluate the performance of the line segmentation algorithm
