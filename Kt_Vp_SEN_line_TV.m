function [Kt,Vp,iter,fres]=Kt_Vp_SEN_line_TV(Kt,Vp,sMaps,U1,kU,opt,Kt_full,Vp_full)
% This is the function to alternatively reconstruct Ktrans and Vp using a
% small number of l-BFGS iterations
% Input: Kt: Initial Ktrans guess
%        Vp: Initial Vp guess
%        sMaps: sensitivity maps
%        U1: under-sampling mask
%        kU: under-sampled k-space
%        opt: other parameters
%Output: Kt,Vp, Kt_inter,Vp_inter: reconstructed and intermediate TK maps
%        fres: objective function value across iterations
%        gradres: gradient value across iterations

% Yi Guo, 06/2014

options.MaxIter=opt.Initer;
options.display = 'off';
options.Method ='lbfgs'; % use l-BFGS method in minFunc
options.optTol=5e-4;
%options.progTol=5e-4;
options.useMex=0;
options.inter=0;
options.numDiff=0;
step = 0.1;
exitflag=0;
iter=1;
fres=[];
omega = 0.05;
beta = 0.5 ;
step = 0.5;
count_Kt = 5;
count_Vp = 5;
% while(exitflag==0)
opt2 = opt;
% lambdaB = [5e-4 0 5e-4 0];
lambdaB = [0 0 0 0];
while(iter <3)
% factor_Kt = 2^iter;
% factor_Vp = 2^iter;

factor_Kt = 1;
factor_Vp = 1;

opt.Vp=reshape(Vp,[opt.size(1) opt.size(2)]);
% opt.lambda1=opt.lambdaA(1)/iter;
% opt.lambda2=opt.lambdaA(2)/iter;
if iter >1
    opt.lambda1=0;
    opt.lambda2=0;
else
    opt.lambda1=opt.lambdaA(1)/(factor_Kt);
    opt.lambda2=opt.lambdaA(2)/(factor_Kt);
end

opt2.lambda1 = lambdaB(1);
if iter==1
    var = find_psnr(Kt,Kt_full)
end
Kt_z = denoised_Kt(Kt) ;
% find_psnr(Kt_z,Kt_full)
% Kt_U = Kt;
 %psnr_val(iter,2) = find_psnr(Kt_z,Kt_full);
% [Kt,f,exitflag,output]=minFunc(@Ktrans2sig_sen_WT,Kt(:),options,sMaps,U1,kU,opt);
cost_val = zeros(20,2);
for count=1:1:count_Kt
    [Kt,cost_val(count,1),cost_val(count,2)]=sd_on_ktrans_search_TV(Kt(:),Kt_z,sMaps,U1,kU,opt,iter,step,beta,omega,opt2);
%     break;
    temp = 2;

end
var_Kt = find_psnr(Kt,Kt_full);

% fres=[fres;output.trace.fval];
%  psnr_val(iter,1) = find_psnr(Kt,Kt_full)
opt.Kt=reshape(Kt,[opt.size(1) opt.size(2)]);
if iter >1
    opt.lambda1=0;
    opt.lambda2=0;
else
    opt.lambda1=opt.lambdaA(1)/(factor_Kt);
    opt.lambda2=opt.lambdaA(2)/(factor_Kt);
end
opt2.lambda1 = lambdaB(3);
% if iter > 5  && iter < 10
%     opt.lambda2 = opt.lambda2/10 ;
% elseif iter > 10
%     opt.lambda2 = opt.lambda2/10;
% end
% [Vp,f,exitflag,output]=minFunc(@Vp2sig_sen_WT,Vp(:),options,sMaps,U1,kU,opt);
cost_val_vp = zeros(20,2);
Vp_z = denoised_Vp(Vp) ;
% find_psnr(Vp,Vp_full)
% find_psnr(Vp_z,Vp_full)

for count=1:1:count_Vp
    [Vp,cost_val_vp(count,1),cost_val_vp(count,2)]=sd_on_vp_search_TV(Vp(:),Vp_z,sMaps,U1,kU,opt,iter,step,beta,omega,opt2);
%     break;
    temp = 2;

end
%  psnr_val(iter,1) = find_psnr(Vp,Vp_full)

var_Vp = find_psnr_Vp(Vp,Vp_full);
var = [var_Kt  var_Vp];
var
% fres=[fres;output.trace.fval];

if opt.plot
    imagesc(real(cat(2,opt.Kt,opt.Vp/2)),[0 0.4]);axis image; axis off; 
    title(['Ktrans, Vp/2, Out iter=',num2str(iter)]); colorbar; 
    drawnow;
end

iter=iter+1;
% iter
end

Kt=reshape(Kt,[opt.size(1) opt.size(2)]);
Vp=reshape(Vp,[opt.size(1) opt.size(2)]);

% var = find_psnr(Kt,Kt_full)

fprintf("slice over")

end


function val = find_psnr(Kt,Kt_r1)

rmse_roi_x = 1: 256;
rmse_roi_y = 1: 256; %slightly bigger



Kt_a = real(Kt);
Kt_b = real(Kt_r1);
Kt_a(Kt_a>0.4) = 0.4 ;
Kt_a(Kt_a<0.0) = 0.0 ;

Kt_b(Kt_b>0.4) = 0.4 ;
Kt_b(Kt_b<0.0) = 0.0 ;

value =psnr(Kt_a(rmse_roi_x,rmse_roi_y),Kt_b(rmse_roi_x,rmse_roi_y));
value2 = ssim(Kt_a(rmse_roi_x,rmse_roi_y),Kt_b(rmse_roi_x,rmse_roi_y));
val = [value, value2];
end

function val = find_psnr_Vp(Vp,Vp_r1)

rmse_roi_x = 1: 256;
rmse_roi_y = 1: 256; %slightly bigger



Vp_a = real(Vp);
Vp_b = real(Vp_r1);
Vp_a(Vp_a>0.8) = 0.8 ;
Vp_a(Vp_a<0.0) = 0.0 ;

Vp_b(Vp_b>0.8) = 0.8 ;
Vp_b(Vp_b<0.0) = 0.0 ;

value =psnr(Vp_a(rmse_roi_x,rmse_roi_y),Vp_b(rmse_roi_x,rmse_roi_y));
value2 = ssim(Vp_a(rmse_roi_x,rmse_roi_y),Vp_b(rmse_roi_x,rmse_roi_y));
val = [value, value2];
end

function Kt_z = denoised_Kt(Kt)
    Kt_abs = real(Kt) ;
%     Kt_abs(Kt_abs<0)= 0;
%     Kt_abs(Kt_abs>0.4)= 0.4;
    save('./network/input_Kt','Kt_abs')
    system('./network/test_Kt_new.sh') ;
    A = load('./network/output_Kt.mat');
    
    Kt_z = squeeze(A.recon);
    temp = 2;

end

function Vp_z = denoised_Vp(Vp)
    Vp_abs = real(Vp) ;
%     Vp_abs(Vp_abs<0)= 0;
%     Vp_abs(Vp_abs>0.8)= 0.8;
    save('./network/input_Vp','Vp_abs')
    system('./network/test_Vp.sh') ;
    A = load('./network/output_Vp.mat');
    
    Vp_z = squeeze(A.recon);
    temp = 2;

end