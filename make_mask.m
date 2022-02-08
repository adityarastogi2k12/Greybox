function mask_total = make_mask(K,R)

[kx,ky,kz,kt,c] = size(K);

if R ==50
    pvalue = 11;
elseif R == 100
    pvalue = 17;
else
    pvalue = 7;
end
mask_total = zeros(kx,ky,kt);
for i = 1:1:kt
    
    res=[kx,ky]; %resolution
    usf=1.0/R; %under sampling factor

    pdf = genPDF(res,pvalue,usf); %density compensation
    mask = genSampling(pdf,2,5);
    mask(1) = 1;
    mask_total(:,:,i) = mask;
end

mask_total = repmat(mask_total,[1,1,1,c]);

mask_total = permute(mask_total,[1,2,5,3,4]);