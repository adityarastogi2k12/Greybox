clc; clear; close all;

mask_total = zeros(100,256,256);
for i = 1:1:100
    nfold=1; 
    res=[256,256]; %resolution
    usf=1.0/nfold; %under sampling factor

    pdf = genPDF(res,7,usf); %density compensation
    mask = genSampling(pdf,2,5);
    mask(1) = 1;
    mask_total(i,:,:) = mask;
end
% imshow(fftshift(mask))

imshow(mask)


%%
% 
% tstMask = mask_total(91:100,:,:);
% trnMask = mask_total(1:90,:,:);
% save('R12','trnMask','tstMask')