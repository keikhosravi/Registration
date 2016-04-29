
clc
clear
close all
tic
he = imread('Y:\Adib\ConklinData\Cropped H&E images for boundary determination\HE\STX_0027 distal normal.tif');
% 
% SHG = imread('Z:\Adib\ConklinData\Cropped H&E images for boundary determination\SHG\STX_0150-01 SHG.tif');
% SHG=im2double(SHG);
%figure;imshow(he)
% max_SHG=max(max(SHG));
% SHG_adj = imadjust(SHG,[0.05 max_SHG],[0 1]);
% SHG_BW=im2bw(SHG_adj,0.001);
% figure;imshow(SHG_BW)
% figure;imshowpair(SHG_adj,he)

S = decorrstretch(he,'tol',0.01);
figure, imshow(he), title('Original Image')
figure,  imshow(S), title('Enhanced Image')

figure;imshow(he), title('H&E image');

cform = makecform('srgb2lab');
lab_he = applycform(he,cform);
ab_gray=im2double(rgb2gray(he));
figure;imshow(ab_gray)
[m n]=size(ab_gray);
  
%  H = fspecial('log',6,0.1);
% deblurred = imfilter(blurred,H,'replicate');
% figure,imshow(deblurred);
H = fspecial('disk',10);
% H = fspecial('gaussian', 10, 100);
    for j=1:3
        k=padarray(S(:,:,j),[70 70],'symmetric');
        k2=histeq(k);
        k1(:,:,j)= imfilter(imfilter(k2,H),H);
        k3(:,:,j)=k1(71:70+m,71:70+n,j);
        
    end
figure, imshow(k3)
%   ab = double(ab_gray);
%   ab = im2double(he(:,:,3));figure; imshow(ab)
    ab=double(k3);
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,nrows*ncols,3);

    nColors =4;
    % repeat the clustering 3 times to avoid local minima
    [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean','Replicates',3);

    pixel_labels = reshape(cluster_idx,nrows,ncols);
    figure;imshow(pixel_labels,[]), title('image labeled by cluster index');

    segmented_images = cell(1,3);
    rgb_label = repmat(pixel_labels,[1 1 3]);
    mean_cluster_intensity=zeros(nColors,1);

        for k = 1:nColors
            color = k3;
            color(rgb_label ~= k) = 0;
            segmented_images{k} = color;
            mean_cluster_intensity(k,1)=mean(nonzeros(rgb2gray(cell2mat(segmented_images(k)))));
            figure;imshow(segmented_images{k})
            tit=['objects in cluster ' num2str(k)];
            title(tit);

        end



     mean_cluster_value = mean(cluster_center,2);
    [tmp, idx] = sort(mean_cluster_value);
    [temp1 idx1]=sort(mean_cluster_intensity);

    cluster_val=zeros(nColors,1);
    for k=1:nColors
        cluster_val(k,1)=find(idx==k)*find(idx1==k);
    end
    [temp2 idx2]=sort(cluster_val);

    blue_cluster_num = idx2(1);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     epith_cell=im2double(cell2mat(segmented_images(blue_cluster_num))+cell2mat(segmented_images(5)));
        epith_cell=im2double(cell2mat(segmented_images(blue_cluster_num)));

% epith_cell=im2double((segmented_images{idx2(2)}));
    figure;imshow(epith_cell)
    
    epith_cell_BW=im2bw(rgb2gray(epith_cell),0.001);
     se = strel('disk',7);
     epith_cell_BW_open = imdilate(epith_cell_BW,se);
     figure;imshowpair(epith_cell_BW_open,he)
     
%      se = strel('disk',10);
%      epith_cell_BW_open = imopen(epith_cell_BW,se);
%      figure;imshowpair(epith_cell_BW_open,he)
     
     BWx= imfill(epith_cell_BW_open,'holes');
     figure;imshowpair(BWx,he)

     BWy=bwareaopen(~BWx,30000);
     BWz=bwareaopen(~BWy,5000);

     
          figure;imshowpair(BWz,he)
closed_mask= imclose(BWz,se);

          figure;imshowpair(closed_mask,he)

     
     
     
     %%%%%%%%%%%color thresholding for finding collagen
     