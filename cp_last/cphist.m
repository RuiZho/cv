function [x_mean,y_mean,C_hist] = cphist(Rec, bw, I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is used to compute the circular projection histogram of
%%% the patterns we found from segmentation.
%%%
%%% Input:  I: the image we hope to segment.
%%%
%%%         bw: the binary value of the image
%%%
%%%         Rec: an k*4 matrix, k is the number of patterns we segmented
%%%         from image I. Each row of Rec contain the position of the
%%%         rectangle that circle the pattern we found.
%%%
%%% Output: C_hist: an 360*k matrix, k is the number of patterns we
%%%         segmented from image I. Each column of C_hist is the circular
%%%         projection histogram of a pattern. It contains 360 values since
%%%         we compute the value every 1 degree.
%%%         x_mean & y_mean: the centroid position of each pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n] = size(I);
figure
for i=1:size(Rec,1)
    % Reconstruct the pattern from Rec.
    im=I(max(floor(Rec(i,2))-1,1):min(floor(Rec(i,2))+Rec(i,4)+1,m),max(floor(Rec(i,1))-1,1):min(floor(Rec(i,1))+Rec(i,3)+1,n),1:3);
    b=bw(max(floor(Rec(i,2))-1,1):min(floor(Rec(i,2))+Rec(i,4)+1,m),max(floor(Rec(i,1))-1,1):min(floor(Rec(i,1))+Rec(i,3)+1,n));
    gim = rgb2gray(im);
    % Delete unwanted part of the pattern, such as small parts from other
    % patterns.
    imLabel = bwlabel(b);
    stats = regionprops(imLabel,'Area');
    area = cat(1,stats.Area);
    index = find(area == max(area));
    b = ismember(imLabel,index);
    gim((b~=1)) = 0;
    % Denoise
    mi = mean(gim(b==1));
    if mi<0.55
        gim(gim>0.65) = 0.1;
    end
    % Roughly compute the shape of the pattern
    gb=edge(b,'canny',0.15);
    subplot(2,size(Rec,1)/2,i);
    imshow(im);
    [x, y] = find(b == 1);
    maxlen = 0;
    minlen = 100000;
    x_mean = mean(x);
    y_mean = mean(y);
    [x, y] = find(gb == 1);
    sx = size(x,1);
    for ii=1:sx
        if sqrt((x(ii)-x_mean)^2+(y(ii)-y_mean)^2)<minlen
            minlen = sqrt((x(ii)-x_mean)^2+(y(ii)-y_mean)^2);
        end
        for jj=ii+1:sx
            if sqrt((x(ii)-x(jj))^2+(y(ii)-y(jj))^2)>maxlen
                maxlen = sqrt((x(ii)-x(jj))^2+(y(ii)-y(jj))^2);
            end
        end
    end
    % Computer the circular projection histogram
    A=zeros(360,1);
    [x, y] = find(b == 1);
    x_mean = mean(x.*gim(b==1));
    y_mean = mean(y.*gim(b==1));
    sx = size(x,1);
    for ii=1:180
        for j=1:sx
            if abs(y(j)-y_mean-tan(ii/180*pi)*(x(j)-x_mean))/sqrt(1+tan(ii/180*pi)*tan(ii/180*pi)) < 4
                if x(j)>x_mean
                    A(ii) = A(ii)+gim(x(j),y(j));
                else
                    A(ii+180) = A(ii+180)+gim(x(j),y(j));
                end
            end
        end
    end
    % Normalize the circular projection histogram, also consider the shape of
    % the pattern.
    A = A/max(A)+mean(mean(gim))+4*minlen/maxlen;
    C_hist(:,i) = A;
end
end

