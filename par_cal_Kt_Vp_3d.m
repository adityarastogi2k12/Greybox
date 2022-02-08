
function [Kt_r1,Vp_r1,Kt,Vp,Kt_U,Vp_U] = par_cal_Kt_Vp_3d(imgF0,lambdaA,R,sigma)


addpath(genpath('Hemant'))
% temp = imgF0(:,:,1,:);
% a = temp(:) ; 
% a = (a - min(a))./(max(a) - min(a));
% imgAll(:,:,:,1) = reshape(a,size(temp));
% 
% imgAll = permute(imgAll,[1,2,4,3]);
imgAll = imgF0/1000;

clear imgF0
%%
opt.d = sqrt(256*256);
% opt.d = 1;
k = fft2(imgAll)./opt.d;
opt.size=size(k);  
[kx,ky,kz,nt]=size(k);
sMaps = ones([kx,ky,kz]) ;
% sMaps=sMaps(:,:,ns,:,:);
sMaps=reshape(sMaps,[kx ky kz 1]);

if ~exist('R1','var')  % use simulated uniform M0 and R1 if none exists
    M0=5*ones(kx,ky,kz,'single'); %use simulated M0, R1
    R1=0.6*ones(kx,ky,kz,'single');
end

imgF= conj(repmat(sMaps,[1 1 1 nt 1])).*ifft2(k).*sqrt(256*256); % get fully-sampled

clear imgF0
clear imgAll
%% set parameters
opt.wname='db4'; % wavelet parameters
opt.worder={[1 2],[1 2],[1 2]};

opt.R1=R1;
opt.M0=M0;
opt.Sb=repmat(imgF(:,:,:,1),[1 1 1 nt]);  %baseline image
opt.alpha=pi*25/180; %flip angle
opt.TR=0.0038;  %TR

delay=5; % delay frames for contrast injection
tpres=4.8/60; % temporal resolution, unit in seconds!
opt.time=[zeros(1,delay),[1:(nt-delay)]*tpres];
opt.plot=0;  % to plot intermediate images during recon

opt.lambdaA=[0.000 0.000 0.000 0.000 0.000 0.000]; % Kt:TV, Wavelet, Vp: TV, wavelet
opt.Initer=10;  %inter interations
opt.Outiter=10; % outer interations

%%
% opt.lambdaA=[0.0000 0.00  0.00000 0.000]; % Kt:TV, Wavelet, Vp: TV, wavelet
opt.lambdaA=lambdaA; % Kt:TV, Wavelet, Vp: TV, wavelet
opt.Initer=10;  %inter interations
opt.Outiter=10; % outer interations
%% calculate fully-sampled Ktrans and Vp
CONCF = sig2conc2((imgF),R1,M0,opt.alpha,opt.TR);
opt.AIF=SAIF_p(opt.time); % get population-averaged AIF
% % opt.AIF = da.aif(1:10:end);
% % opt.AIF = AIF(floor(linspace(1,length(AIF),nt)),1)'./1 ;
% % opt.AIF(4:end) = opt.AIF(3:end-1);
% % opt.AIF(3) = 0;
% % opt.AIF(opt.AIF==0) = 0.001;
% % opt.AIF = Cp0.*(1-0.4);
[Kt,Vp]=conc2Ktrans(CONCF,opt.time,opt.AIF);
%%
% R = 100; %under-sampling rate

U1 = make_mask(k,R) ;
% numel(U1)./sum(U1(:))
U1(:,:,:,1,:)=1;  % first frame is fully-sampled
% U1(:,:,:,1,:)=1;  % first frame is fully-sampled
kU=k.*U1;

% numel(U1)./sum(U1(:))
if sigma ~=0
    noise = randn(size(kU)) + 1i.*randn(size(kU));
    noise = noise.*sigma./sqrt(2);
    noise(:,:,:,1,:) = 0;
    kU = kU+noise;
end

imgU= conj(repmat(sMaps,[1 1 1 nt])).*ifft2(kU).*opt.d;
CONC1 = sig2conc2(real(imgU),R1,M0,opt.alpha,opt.TR);
[Kt_U,Vp_U]=conc2Ktrans(CONC1,opt.time,opt.AIF);
%%
clear('imgF','imgU','k','M0','R1','U11','CONC1','CONCF','a')
% % % clc
Kt_2=zeros(opt.size(1),opt.size(2),opt.size(3)); 
% Vp_1=Kt_1;
Kt_1 = Kt_U;
Vp_1 = Vp_U ;
Vp_2 = Kt_2;
opt2 = opt;
opt2.lambdaA= [0 0 0 0];
% [Kt_L2,Vp_L2,~,~]=Kt_Vp_SEN_org(Kt_2,Vp_2,sMaps,U1,kU,opt2,Kt);
% % % % tic,
[Kt_1,Vp_1,~,~]=Kt_Vp_SEN_line_TV(Kt_U,Vp_U,sMaps,U1,kU,opt,Kt,Vp);
% % % toc,
Kt_r1=real(Kt_1);  % get real value after recon
Vp_r1=real(Vp_1);
fprintf('End')
