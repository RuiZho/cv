%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is for a single detect. It shows the segmentation and
%%% matching results and it will not store the result. You will have to
%%% manually change the value of title1 and title before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
% Read two images
title1 = '.\cards\001.png';
title2 = '.\cards\002.png';
I1= double(imread(title1))/255;
I2= double(imread(title2))/255;
% Segmentation and find feature discripter
[Rec1, bw1] = segmentation(I1);
[Rec2, bw2] = segmentation(I2);
[xm,ym,B] = cphist(Rec1,bw1,I1);
[~,~,C] = cphist(Rec2,bw2,I2);
BB = colorhist(Rec1,bw1,I1);
CC = colorhist(Rec2,bw2,I2);
% Matching by color histogram
count1 = size(B,2);
count2 = size(C,2);
no = 100000*ones(count1,count2);
for i=1:count1
    for j=1:count2
        c=C(:,j);
        b=B(:,i);
        siz = size(b,1);
        for k=1:siz
            b=[b(2:end);b(1)];
            no(i,j) = min(no(i,j),norm(c-b,1));
        end
    end
end
ns = sort(no(:));
for i=1:size(ns,1)
    [mm(i,1), mm(i,2)] = find(no==ns(i));
end
mm = mm(1:7,:);
% Matching by circular projection histogram
s1=size(BB,3);
s2=size(CC,3);
c=ones(s1,s2)*100;
thr = 0.1;
for i=1:s1
    for j=1:s2
        f=abs(BB(:,:,i)-CC(:,:,j));
        c1(i,j)=sum(sum(f>thr));
        c2(i,j)=max(max(f));
        if(sum(sum(f>thr))<9)
            c(i,j)=norm(f,1);
        elseif(sum(sum(f>thr))<11)
            c(i,j)=norm(f,1)+9;
        elseif(sum(sum(f>thr))<13)
            c(i,j)=norm(f,1)+11;
        end
    end
end
cs = sort(c(:));
for i=1:size(cs,1)
    if i<size(cs,1) && cs(i)==cs(i+1)
        break;
    end
    [m(i,1),m(i,2)] = find(c==cs(i));
end
% Make trade off between two matching results
if c1(m(1,1),m(1,2))>c1(m(2,1),m(2,2))
    tmp = m(1,:);
    m(1,:) = m(2,:);
    m(2,:) = tmp;
end
x=0;
y=0;
for i=1:size(m,1)
    for j=1:7
        if mm(j,1)==m(i,1) && mm(j,2)==m(i,2)
            x=m(i,1);
            y=m(i,2);
            break;
        end
    end
    if x~=0
        break;
    end
end
% Show and store the result
if x~=0
    [mx1,nx1,~] = size(I1);
    [mx2,nx2,~] = size(I2);
    im1=I1(max(floor(Rec1(x,2))-1,1):min(floor(Rec1(x,2))+Rec1(x,4)+1,mx1),max(floor(Rec1(x,1))-1,1):min(floor(Rec1(x,1))+Rec1(x,3)+1,nx1),1:3);
    im2=I2(max(floor(Rec2(y,2))-1,1):min(floor(Rec2(y,2))+Rec2(y,4)+1,mx2),max(floor(Rec2(y,1))-1,1):min(floor(Rec2(y,1))+Rec2(y,3)+1,nx2),1:3);
    figure
    imshow(im1);
    figure
    imshow(im2);
else
    disp('Fail to match');
end