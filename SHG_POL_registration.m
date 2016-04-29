clc
clear
close all

SHG_list = dir('Z:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\SHG1');
POL_path = 'Z:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\POL1\';
SHG_path='Z:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\SHG1\';
%width_mean=zeros(length(dirlist),1);
for fileind = 1:length(SHG_list)
       tic 
%        fileind=18;
    SHG_file_spec=SHG_list(fileind);
    [SHG_pathstr,SHG_name,SHG_ext] = fileparts(SHG_file_spec.name);
    
    if(strcmp(SHG_ext,'.tif'))
       
            SHG_image=[SHG_path SHG_name SHG_ext];
            POL_image=[POL_path SHG_name(1:end-3) 'POL' SHG_ext];
        
    
POL = imread(POL_image);
POL1=fliplr(im2double(POL));
POL=imadjust(POL1,[0.15,0.9],[0,1]);

SHG = imread(SHG_image);
fixed=imadjust(SHG);

% figure;imshow(fixed)
moving=imresize(POL,size(fixed));
% figure;imshow(moving)

[optimizer,metric] = imregconfig('multimodal');
% movingRegisteredDefault = imregister(moving, fixed, 'affine', optimizer, metric);
% figure, imshowpair(movingRegisteredDefault, fixed)
% title('A: Default registration')


disp(optimizer)
disp(metric)
optimizer.InitialRadius = optimizer.InitialRadius/3.5;
movingRegisteredAdjustedInitialRadius = imregister(moving, fixed, 'affine', optimizer, metric);
% figure, imshowpair(movingRegisteredAdjustedInitialRadius, fixed)
% title('Adjusted InitialRadius')
optimizer.MaximumIterations = 700;
movingRegisteredAdjustedInitialRadius700 = imregister(moving, fixed, 'affine', optimizer, metric);
% figure, imshowpair(movingRegisteredAdjustedInitialRadius300, fixed)
% title('B: Adjusted InitialRadius, MaximumIterations = 700, Adjusted InitialRadius.')
tformSimilarity = imregtform(moving,fixed,'similarity',optimizer,metric);
Rfixed = imref2d(size(fixed));
movingRegisteredRigid = imwarp(moving,tformSimilarity,'OutputView',Rfixed);
figure, imshowpair(movingRegisteredRigid, fixed);
title('C: Registration based on similarity transformation model.');
tformSimilarity.T
[movingRegisteredAffineWithIC treg tform]= imreg_new3(moving,fixed,'affine',optimizer,metric,...
'InitialTransformation',tformSimilarity);
% figure, imshowpair(movingRegisteredAffineWithIC,fixed);
% title('D: Registration from affine model based on similarity initial condition.');
% 
% tform= imregtform(moving,fixed,'affine',optimizer,metric,...
% 'InitialTransformation',tformSimilarity);
% movingRegisteredAffineWithIC1=imwarp(moving,tform);
% figure, imshowpair(movingRegisteredAffineWithIC1,fixed);
% title('D: Registration from affine model based on similarity initial condition1.');
%%%%%%%
Rmoving=imref2d(size(moving));
POL_registered=imresize(POL1,size(fixed));
% tform = imregtform(moving, movingRegisteredAffineWithIC, 'affine', optimizer, metric);
B = imwarp(POL_registered,Rmoving,tform,'OutputView',Rfixed);
figure
imshowpair(B, fixed)
% title('registered');

Registered_image=['Z:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\Registered_POL\' SHG_name(1:end-3) 'POL_registered' SHG_ext]; 
imwrite(B,Registered_image);

toc
    end
end

