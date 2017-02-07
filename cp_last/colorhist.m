function Color = colorhist(Rec, bw, I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is used to compute the color histogram of the patterns we
%%% found from segmentation.
%%%
%%% Input:  I: the image we hope to segment.
%%%
%%%         bw: the binary value of the image
%%%
%%%         Rec: an k*4 matrix, k is the number of patterns we segmented
%%%         from image I. Each row of Rec contain the position of the
%%%         rectangle that circle the pattern we found.
%%%
%%% Output: Color: an 16*3*k matrix, k is the number of patterns we
%%%         segmented from image I. For each pattern, we have 16 values for
%%%         each of the R, G and B layer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n] = size(I);
for i=1:size(Rec,1)
    % Reconstruct the pattern from Rec and resize it.
    im=I(max(floor(Rec(i,2))-1,1):min(floor(Rec(i,2))+Rec(i,4)+1,m),max(floor(Rec(i,1))-1,1):min(floor(Rec(i,1))+Rec(i,3)+1,n),1:3);
    b=bw(max(floor(Rec(i,2))-1,1):min(floor(Rec(i,2))+Rec(i,4)+1,m),max(floor(Rec(i,1))-1,1):min(floor(Rec(i,1))+Rec(i,3)+1,n));
    im=imresize(im,[250,250],'nearest');
    b=imresize(b,[250,250],'nearest');
    gim = rgb2gray(im);
    % Delete unwanted part of the pattern, such as small parts from other
    % patterns.
    [mm,nn] = size(b);
    imLabel = bwlabel(b);
    stats = regionprops(imLabel,'Area');
    area = cat(1,stats.Area);
    index = find(area == max(area));
    b = ismember(imLabel,index);
    gim((b~=1)) = 1;
    gim(gim>0.75)=1;
    gim(gim<0.1)=0;
    % Denoise
    for ii=1:mm
        for j=1:nn
            if (abs(gim(ii,j)-1)<0.2 )
                im(ii,j,:)=[1 1 1];
            end
        end
    end
    % Compute color histogram
    [mm, nn, ~]=size(im);
    countR=1;countG=1;countB=1;
    Red=0;Green=0;Blue=0;
    for ii=1:mm
        for jj=1:nn
            if(im(ii,jj,1)~=1)
                Red(countR,1)=im(ii,jj,1);
                countR=countR+1;
            end
            if(im(ii,jj,2)~=1)
                Green(countG,1)=im(ii,jj,2);
                countG=countG+1;
            end
            if(im(ii,jj,3)~=1)
                Blue(countB,1)=im(ii,jj,3);
                countB=countB+1;
            end
        end
    end
    [y(:,1),~] = imhist(Red,16);
    [y(:,2),~] = imhist(Green, 16);
    [y(:,3),~] = imhist(Blue, 16);
    y(:,1)=y(:,1)./size(Red,1);
    y(:,2)=y(:,2)./size(Green,1);
    y(:,3)=y(:,3)./size(Blue,1);
    Color(:,:,i)=y;
end

end

