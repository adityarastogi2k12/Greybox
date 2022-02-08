
close all; 
clc;
clear ;
m = [1,5,6,8,10];
V = 1;
addpath(genpath('minFunc_2012'));
addpath(genpath('ComputingFunction'));
addpath(genpath('Training_data'));
%%

data_name  = 'braindataPATC' ;
load(data_name) ;

zdim = [1:14];
[nx,ny,nz,nt] = size(imgs2(:,:,zdim,:));
imgF0 = imgs2(:,:,zdim,:);
clear('imgs')
% imgF0 = parallel.pool.Constant(imgF0);
%%
Kt =zeros(nx,ny,nz);
Vp = Kt ;
Kt_U = Kt ;
Vp_U = Kt ;
Kt_R= Kt ;
Vp_R= Kt ;
Kt_L2 = Kt;
Vp_L2 = Kt;
Kt_NN = Kt;
Vp_NN = Kt; 
Kt_TV = Kt;
Vp_TV = Vp ;

R = 20;
tic
lambdaA=[0.0 0.00  0.0  00]; % Kt:TV, Wavelet, Vp: TV, wavelet
% a= parallel.pool.Constant(lambdaA);
% b = parallel.pool.Constant(R);
a = lambdaA ;
a2 = [0 0 0 0];
b = R;
sigma = 0;
%%
% lambdaA=[0.1 0.00 0.01 00]; % Kt:TV, Wavelet, Vp: TV, wavelet
lambdaA=[5 0.00 5 00]; % Kt:TV, Wavelet, Vp: TV, wavelet
% best results so farfor pat 35 with lambdaA=[0.1 0.00  0.1  00] count_Kt = 5;
% count_Vp = 5; lambdaB = [5e-3 0 5e-4 0]; iter <11  factor_Kt = iter; factor_Vp = iter;

a = lambdaA ;
% for i = 1:1:length(zdim)
for i = 11
    i
    [Kt_R(:,:,i),Vp_R(:,:,i),Kt(:,:,i),Vp(:,:,i),Kt_U(:,:,i),Vp_U(:,:,i)] = par_cal_Kt_Vp_3d(imgF0(:,:,i,:),a,b,sigma);
    
end
toc

%%
% save('PatC_R20.mat','Kt','Kt_R','Kt_U','Vp','Vp_R','Vp_U')
