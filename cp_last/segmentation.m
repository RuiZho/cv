function [Rec, bw] = segmentation(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is used to segment image and extract the pattern in the
%%% iamge.
%%%
%%% Input:  I: the image we hope to segment.
%%%
%%% Output: Rec: an k*4 matrix, k is the number of patterns we segmented
%%%         from image I. Each row of Rec contain the position of the
%%%         rectangle that circle the pattern we found.
%%%
%%%         bw: the binary value image of I.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I1=edge(rgb2gray(I),'canny',0.15);
[m, n]=size(I1);
% Find the closed region in the image.
bw = imclose(I1, strel('disk', 2));
bw = imfill(bw, 'holes');
bw = imopen(bw, strel('disk', 2));
a=regionprops(bw,'Area');
e=edge(bw,'Canny');
% If one of the regions is too big, abandon it.
if (max([a.Area])/(m*n)>0.5)
    I1(e == 1) = 0;
end
% Find the closed regions in the image.
bw = imclose(I1, strel('disk', 2));
bw = imfill(bw, 'holes');
bw = imopen(bw, strel('disk', 2));
bw_img=bw;
img_reg = regionprops(bw_img,  'area', 'boundingbox','image');
% Find the boundaries of the regions.
areas = [img_reg.Area];
rects = cat(1,img_reg.BoundingBox);
% Find the position of the patterns, if a pattern is too small, abandon it.
count=1;
for i=1:size(areas,2)
    if(areas(1,i)/(m*n)>0.0007)
        Rec(count,:)=rects(i,:);
        count=count+1;
    end
end
end
