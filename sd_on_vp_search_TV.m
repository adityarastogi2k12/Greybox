function [Vp,cost1,cost2]=sd_on_vp_search_TV(Vp,Vp_z,sMaps,U1,kU,opt,iter,step,beta,omega,opt2)
% Least-square equation and L1 wavelet and TV with direvative calculation
% from Vp to Kspace

% Yi Guo, 06/12/2014

%% Calculate cost function
Rcs=3.8;
alpha=opt.alpha;
TR=opt.TR;

% alpha = 0.1;
% beta = 0.1 ;
% step = 0.1;
Vp=reshape(Vp,[opt.size(1) opt.size(2)]);
[cost1,CAIF,CR,S] = get_cost(Vp,kU,opt,U1,alpha,sMaps) ;
if opt2.lambda1 ~= 0
    cost2 = get_cost_TV(Vp,opt2);
    cost1 = cost1+ cost2;
end
cost1 = cost1 + opt.lambda1.*norm((Vp-Vp_z),2) ;
g = get_grad(S,kU,CAIF,opt,CR,Rcs,TR,sMaps,alpha) ;
g1 = g + opt.lambda1.*(Vp-Vp_z) ;
% sprintf('g1_norm= %0.5d, g2_norm = %0.5d',norm(g1(:)),norm(opt.lambda1.*(Vp-Vp_z)))
if opt2.lambda1 ~= 0
    g2 = get_grad_TV(Vp,opt2);
    g1 = g1+ g2;
end

Vp_temp = Vp - step.*g1;
[cost1_temp,~,~,~] = get_cost(Vp_temp,kU,opt,U1,alpha,sMaps) ;
Delta=  omega*(norm(g1(:)))^2;
if opt2.lambda1 ~= 0
    cost2_temp = get_cost_TV(Vp_temp,opt2);
    cost1_temp = cost1_temp+ cost2_temp;
end
cost1_temp = cost1_temp + opt.lambda1.*norm((Vp_temp-Vp_z),2) ;
while cost1_temp > cost1 - step*Delta
    step = step*beta;
    if step < 1e-5
        step = 1e-5;
        Vp_temp = Vp - step.*g1;
        break
    end
    Vp_temp = Vp - step.*g1;
    [cost1_temp,~,~,~] = get_cost(Vp_temp,kU,opt,U1,alpha,sMaps) ;
    cost1_temp = cost1_temp + opt.lambda1.*norm((Vp_temp-Vp_z),2) ;
    if opt2.lambda1 ~= 0
        cost2_temp = get_cost_TV(Vp_temp,opt2);
        cost1_temp = cost1_temp+ cost2_temp;
    
    end
end
Vp = Vp_temp;
cost2 = 0;
% if opt.lambda1 ~= 0
%     [cost_reg,g_reg] = get_grad_reg(Vp)   ;
%     cost2 = opt.lambda1.*cost_reg;
%     g2 = opt.lambda1.*g_reg;
%     sprintf('data cost =%0.5d and reg cost = %0.5d, g1_max= %0.5d, g2_max = %0.5d',cost1,cost2,norm(g1(:)),norm(g2(:)))
%     Vp = Vp -step.*g2 ;
% else
%     cost2 = 0;
% end

% 
% if opt.lambda1~=0 % calculate lambda1*||TVx||1
% % TV=compTx2d(Vp,opt);
% % cost2=opt.lambda1*sum(abs(TV(:)));
% Vp_abs = real(Vp) ;
% save('../input_Vp','Vp_abs')
% [cost_reg,grad_reg] = get_grad_reg(Vp_abs) ;
% cost2 = cost_reg ;
% else
%     cost2=0;
% end

% if opt.lambda2~=0  % calculate lambda2*||Wx||1
% [wc,fsize]=compWx(Vp,opt);
% cost3=opt.lambda2*sum(abs(wc(:)));
% else
%     cost3=0;
% end

% if opt.lambda2~=0  % calculate lambda2*||Wx||1
%     % [wc,fsize]=compWx(Vp,opt);
% %     wc = dct2(Vp);
%     wc  = Vp;
%     cost3=opt.lambda2*sum(abs(wc(:)));
% else
%     cost3=0;
% end
% 
% 
% cost=cost1+cost2+cost3;

%% calculate gradient of cost function


% if opt.lambda1~=0
% %     TV=compTx2d(Vp,opt);
% % u=1e-8;
% % W=sqrt(conj(TV).*TV+u);
% % W=1./W;
% % g2=opt.lambda1.*compThx2d(W.*TV,opt);
% g2 = opt.lambda1.*grad_reg ;
% else
%     g2=0;
% end
% 
% % if opt.lambda2~=0
% % %[wc,fsize]=compWx(Vp,opt);
% % u=1e-8;
% % W1=sqrt(conj(wc).*wc+u);
% % W1=1./W1;
% % g3=opt.lambda2.*compWhx(W1.*wc,opt,fsize);
% % else
% %     g3=0;
% % end
% 
% if opt.lambda2~=0
% %[wc,fsize]=compWx(Kt,opt);
%     u=1e-8;
%     W1=sqrt(conj(wc).*wc+u);
%     W1=1./W1;
%     % g3=opt.lambda2.*compWhx(W1.*wc,opt,fsize);
%     g3 = reshape(W1.*wc,[opt.size(1) opt.size(2)]);
% %     g3 = opt.lambda2.*idct2(g3);
%     g3 = opt.lambda2.*g3;
% else
%     g3=0;
% end
%  sprintf('g1_norm= %0.5d, g2_norm = %0.5d',norm(g1(:)),norm(opt.lambda1.*(Vp-Vp_z)))
% 
% step_size = 0.001;
% if (iter > 50) && (iter <100)
%     step_size = step_size/2 ;
% elseif iter > 100
%     step_size = step_size/4 ;
% end
% 
% Vp_1 = Vp - step_size.*g1;
% Vp = Vp_1 -step_size.*g2 ;
% % imagesc(abs(Vp),[0 0.4])
% Vp = Vp ;
% 
% grad=g1+g2+g3;
% 
% grad=(grad(:));

end

function [cost,CAIF,CR,S] = get_cost(Vp,kU,opt,U1,alpha,sMaps)
    nt=opt.size(4);
    Vp=reshape(Vp,[opt.size(1) opt.size(2)]);
    nt=opt.size(4);

    [CR,~,CAIF]=Ktrans2conc(opt.Kt,Vp,opt.time,opt.AIF);
    S = conc2sig(CR,opt.R1,opt.M0,opt.Sb,alpha,opt.TR);
    S=U1.*fft2(repmat(sMaps,[1 1 1 nt]).*repmat(S,[1 1 1 1]))./opt.d;
    cost=0.5*sum(abs(S(:)-kU(:)).^2);
end

function cost2 = get_cost_TV(Vp,opt)


    TV=compTx2d(Vp,opt);
    cost2=opt.lambda1*sum(abs(TV(:)));

end

function [g1] = get_grad(S,kU,CAIF,opt,CR,Rcs,TR,sMaps,alpha)
    nt=opt.size(4);
    Rt=CR*Rcs+repmat(opt.R1,[1 1 1 nt]);

    E1=exp(-Rt*TR);
    dBdC=(TR*E1.*(1-cos(alpha).*E1)-(1-E1).*cos(alpha).*TR.*E1)./((1-cos(alpha).*E1).^2);
    dBdC=dBdC.*Rcs;

    g=dBdC.*repmat(opt.M0,[1 1 1 nt]).*sin(alpha).*(repmat(conj(sMaps),[1 1 1 nt]).*ifft2(S-kU).*opt.d);

    g1=zeros(opt.size(1),opt.size(2));

    for it=1:nt
    g1=g1+g(:,:,:,it).*CAIF(it);
    end

end

function [g2] = get_grad_TV(Vp,opt)
    u=1e-8;
    TV=compTx2d(Vp,opt);
    W=sqrt(conj(TV).*TV+u);
    W=1./W;
    g2=opt.lambda1.*compThx2d(W.*TV,opt);

end

function [cost,grad_reg] = get_grad_reg(Vp)
    Vp_abs = real(Vp) ;
    save('../input_Vp','Vp_abs')
    system('../test_Vp.sh') ;
    A = load('../grad_Vp.mat');
    cost = A.cost;
    grad_reg = squeeze(A.grad);
    temp = 2;

end
