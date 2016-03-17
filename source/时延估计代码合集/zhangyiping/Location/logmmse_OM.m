function parameters=logmmse_OM(ns_ps,parameters)
len = parameters.len;
Gmin=parameters.Gmin;  % 无音时假设非零但很小，Gmin为小值
noise_ps=parameters.noise_ps;
eps_min=parameters.eps_min;
beta=parameters.beta;  %式3-65时域平滑参数
psi=parameters.psi;
psi_min=parameters.psi_min;  %式3-67、3-68限值
psi_max=parameters.psi_max;
psi_pmin=parameters.psi_pmin;
psi_pmax=parameters.psi_pmax; 
psi_fr_old=parameters.psi_fr;
psi_peak=parameters.psi_peak;
qmax=parameters.qmax;
alpha=parameters.alpha; %式3-63平滑因子
GH1=parameters.GH1;
gamma_old = parameters.gamma; %前一帧后验信噪比
eps_old = parameters.eps;

Sf = zeros(len,1);
Sf(1)=ns_ps(1);   %带噪声能量谱的频域平滑
Sf(end)=ns_ps(end);
for f=2:len-1
    Sf(f)=0.25*ns_ps(f-1)+0.5* ns_ps(f)+0.25* ns_ps(f+1);   %带噪能量谱使用3阶hanning窗进行平滑
end
% sig=sqrt(ns_ps);  %带噪信号频谱
gamma=min(Sf./noise_ps,40);  % 后验信噪比，限值后验性噪比最大值为16dB？？
% eps=alpha*X2_prev./noise_mu2 + (1-alpha)*max(gamma-1,0);   %先验信噪比，式3-63   
eps=alpha*(GH1.^2).*gamma_old + (1-alpha)*max(gamma-1,0);
eps=max(eps_min,eps);  % 控制先验信噪比大于 -25 dB？？

A=eps./(1+eps);
v=A.*gamma;
ei_v=0.5*expint(v);
GH1=A.*exp(ei_v);    %计算GH1，可使用噪声估计结果？？？
% GH1=max(1,GH1); 

% ---估算条件无音概率 P(Ho) ---------------

len2=len/2+1; 
psi=beta*psi+(1-beta)*eps_old(1:len2); %公式3-65先验信噪比平滑       
psi_local=smoothing(psi,1);  %加窗平滑，先验信噪比的局部和全局平滑值
psi_global=smoothing(psi,5);  
C=log10(psi_max/psi_min); 

Plocal=zeros(len2,1);   % 计算 P_local
imax=find(psi_local>=psi_max);
Plocal(imax)=1;
ibet=find(psi_local>psi_min & psi_local<psi_max);
Plocal(ibet)=log10(psi_local(ibet)/psi_min)/C;
    
    
Pglob=zeros(len2,1);   % 计算 P_global
imax=find(psi_global>=psi_max);
Pglob(imax)=1;
ibet=find(psi_global>psi_min & psi_global<psi_max);
Pglob(ibet)=log10(psi_global(ibet)/psi_min)/C;
    
psi_fr=mean(psi);  % 计算 Pframe
if  psi_fr>psi_min
    if psi_fr>psi_fr_old
       Pframe=1;
       psi_peak=min(max(psi_fr,psi_pmin),psi_pmax);
    else
        if psi_fr <=psi_peak*psi_min
           Pframe=0;
        elseif psi_fr>= psi_peak*psi_max
            Pframe=1;
        else
            Pframe=log10(psi_fr/psi_peak/psi_min)/C;
        end
     end
else
    Pframe=0;
end

q2 = 1- Plocal.*Pglob*Pframe;  % 先验无音概率
% q2= min(qmax,q2); %限值无音概率低于qmax<1
q = [q2; flipud(q2(2:len2-1))]; %复制并反向，形成FFT分析长度 

pSAP = (1-q)./(1-q+q.*(1+eps).*exp(-v)); % P(H1 | Yk)
iqmax=find(q>qmax);
pSAP(iqmax)=0;


if pSAP > 0.8
   Gmin2 = Gmin.^(1-pSAP)*2;
else
   Gmin2 = Gmin.^(1-pSAP); 
end
Gcohen=(GH1.^pSAP).*Gmin2;
% sig = sig.*Gcohen; 
X2=Gcohen.^2;

parameters.X2=X2;
parameters.psi=psi;
parameters.psi_fr=psi_fr; 
parameters.psi_peak=psi_peak;    
parameters.sigout=Gcohen;  
parameters.GH1=GH1; 
parameters.gamma = gamma;
parameters.eps = eps;

   
    


 



