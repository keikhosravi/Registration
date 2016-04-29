clc
clear
close all

SHG_list = dir('Y:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\SHG');
SHG_path='Y:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\SHG\';

for fileind = 1:length(SHG_list)
       tic 
    SHG_file_spec=SHG_list(fileind);
    [SHG_pathstr,SHG_name,SHG_ext] = fileparts(SHG_file_spec.name);
    
    if(strcmp(SHG_ext,'.tif'))
        
     SHG_image=[SHG_path SHG_name SHG_ext];
     B=imread(SHG_image);
       
     SHG_name1=['Core' SHG_name(1:end-5) '_ROI' SHG_name(end-3) '_SHG'];
     
     Registered_image=['Y:\Adib\LOCI-LC-PolScope\pancreatic TMA\Registered\SHG\' SHG_name1 SHG_ext]; 
imwrite(B,Registered_image);

    end
end
