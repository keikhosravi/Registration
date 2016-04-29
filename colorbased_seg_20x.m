
    clc
    clear
    close all
    tic
    he = imread('Z:\Cole\Teresa TMA SHG ROIs\Cropped HE\1C HE_registered_ROI2.tif');

    figure;imshow(he), title('H&E image');

    cform = makecform('srgb2lab');
    lab_he = applycform(he,cform);
    ab_gray=im2double(rgb2gray(he));
    figure;imshow(ab_gray)
  
%  H = fspecial('log',6,0.1);
% deblurred = imfilter(blurred,H,'replicate');
% figure,imshow(deblurred);
H = fspecial('disk',20);
% H = fspecial('gaussian', 10, 100);
    for j=1:3
        k=he(:,:,j);
        k2=histeq(k);
        k1(:,:,j)= imfilter(imfilter(k2,H),H);
        
    end
figure, imshow(k1)
%   ab = double(ab_gray);
%   ab = im2double(he(:,:,3));figure; imshow(ab)
    ab=double(k1);
    nrows = size(ab,1);
    ncols = size(ab,2);
    ab = reshape(ab,nrows*ncols,3);

    nColors =5;
    % repeat the clustering 3 times to avoid local minima
    [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean','Replicates',3);

    pixel_labels = reshape(cluster_idx,nrows,ncols);
    figure;imshow(pixel_labels,[]), title('image labeled by cluster index');

    segmented_images = cell(1,3);
    rgb_label = repmat(pixel_labels,[1 1 3]);
    mean_cluster_intensity=zeros(nColors,1);

        for k = 1:nColors
            color = k1;
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
    epith_cell=im2double(cell2mat(segmented_images(blue_cluster_num))+cell2mat(segmented_images(idx2(2))));
% epith_cell=im2double((segmented_images{idx2(2)}));
    figure;imshow(epith_cell)
    
    epith_cell_BW=im2bw(rgb2gray(epith_cell),0.001);
     se = strel('disk',15);
     epith_cell_BW_open = imopen(epith_cell_BW,se);
     figure;imshowpair(epith_cell_BW_open,he)


    
    % L = lab_he(:,:,1);
    % blue_idx = find(pixel_labels == blue_cluster_num);
    % L_blue = L(blue_idx);
    % is_light_blue = im2bw(L_blue,graythresh(L_blue));
    %
    % nuclei_labels = repmat(uint8(0),[nrows ncols]);
    % nuclei_labels(blue_idx(is_light_blue==false)) = 1;
    % nuclei_labels = repmat(nuclei_labels,[1 1 3]);
    % blue_nuclei = he;
    % blue_nuclei(nuclei_labels ~= 1) = 0;
    % figure;imshow(blue_nuclei), title('blue nuclei');
    
    Collagen_red=cell2mat(segmented_images(idx2((end-1))));
    Collagen_BW=im2bw(im2double(rgb2gray(Collagen_red)),0.001);
    figure;imshow(Collagen_red)
    se = strel('disk',4);
    Collagen_BW_dilate = imopen(Collagen_BW,se);
    figure;imshow(Collagen_BW_dilate)
     se = strel('disk',20);
    Collagen_BW_close = imclose(Collagen_BW_dilate,se);
    figure;imshow(Collagen_BW_close)
    
    epith_cell=im2double(cell2mat(segmented_images(blue_cluster_num))+cell2mat(segmented_images(idx2(2))));
% epith_cell=im2double((segmented_images{idx2(2)}));
    figure;imshow(epith_cell)
        figure;imshowpair(epith_cell,Collagen_BW_close)

    macro_cell=im2double(cell2mat(segmented_images(blue_cluster_num)));
%         figure;imshow(macro_cell)

    ab1=double(macro_cell);
    nrows = size(ab1,1);
    ncols = size(ab1,2);
    ab1 = reshape(ab1,nrows*ncols,3);

    nColors = 4;
    % repeat the clustering 3 times to avoid local minima
    [cluster_idx1, cluster_center1] = kmeans(ab1,nColors,'distance','sqEuclidean','Replicates',3);

    pixel_labels = reshape(cluster_idx1,nrows,ncols);
    figure;imshow(pixel_labels,[]), title('image labeled by cluster index');

    segmented_images = cell(1,3);
    rgb_label = repmat(pixel_labels,[1 1 3]);
    mean_cluster_intensity=zeros(nColors,1);

        for k = 1:nColors
            color = macro_cell;
            color(rgb_label ~= k) = 0;
            segmented_images{k} = color;
            mean_cluster_intensity(k,1)=mean(nonzeros(rgb2gray(cell2mat(segmented_images(k)))));
            figure;imshow(segmented_images{k})
            tit=['objects in cluster ' num2str(k)];
            title(tit);

        end
 mean_cluster_value = mean(cluster_center1,2);
    [tmp, idx] = sort(mean_cluster_value);
    
%     cluster_val=zeros(nColors,1);
%     for k=1:nColors
%         cluster_val(k,1)=find(idx==k)*find(idx1==k);
%     end
%     [temp2 idx2]=sort(cluster_val);

    blue_cluster_num = idx(2);
    
    BW_macro=im2bw(im2double(rgb2gray(segmented_images{blue_cluster_num})),0.001);
    figure;imshow(BW_macro)
    
    BW_macro_disc=bwareaopen(BW_macro,10);
    figure;imshow(BW_macro_disc)
    se=strel('disk',4);
    BW_macro_open=imopen(BW_macro_disc,se);
    figure;imshow(BW_macro_open)
    se=strel('disk',10);
    BW_macro_dil=imdilate(BW_macro_open,se);
    figure;imshow(BW_macro_dil)
    for i=1:3
epith_cell_discard(:,:,i)=im2double(~BW_macro_dil).*epith_cell(:,:,i);
    end
    figure;imshow(epith_cell_discard)

    
    I_gray=im2double(rgb2gray(epith_cell_discard));
    BW=im2bw(I_gray,0.001);
    figure;imshow(BW)
    se=strel('disk',3);
    BW_dil=imdilate(BW,se);
    figure;imshow(BW_dil)
    BW_filled = imfill(BW_dil,'holes');
    figure;imshow(BW_filled)
    BW_disc=bwareaopen(BW_filled,3000);
    figure;imshow(BW_disc)

    for i=1:3
        epith_cell(:,:,i)= im2double(BW_disc).*im2double(he(:,:,i));
    end
    
    figure;imshow(epith_cell)

       
    se = strel('disk',2);
    BW=imerode(BW_disc,se);
    figure;imshowpair(BW_disc,he)
    BW2 = BW.*(~Collagen_BW_dilate);
    figure
    imshow(BW2)
    I_gray_discard=I_gray.*BW2;
    figure
    imshow(I_gray_discard)



    I_filt=imfilter(I_gray_discard,fspecial('average',10));
    % SE = strel('disk', 5, 0);
    %  I_filt = imdilate(I_gray_discard,SE);
    I_filt=imfilter(I_filt,fspecial('gaussian',20,100));
    % I_filt=imfilter(I_filt,fspecial('average',20));

    I_filt = imadjust(I_filt);
    I_filt_BW=im2bw(I_filt,0.01);
    figure; imshow(I_filt)
    figure;imshow(I_filt_BW)


    IM2 = imcomplement(I_filt_BW);
    IM2 = bwareaopen(IM2, 6000);
    IM3 = imcomplement(IM2);
    figure;imshow(IM3)
    figure;imshowpair(IM3,he)


    
    
    toc

