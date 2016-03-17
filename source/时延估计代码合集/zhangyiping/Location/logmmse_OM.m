function parameters=logmmse_OM(ns_ps,parameters)
len = parameters.len;
Gmin=parameters.Gmin;  % ����ʱ������㵫��С��GminΪСֵ
noise_ps=parameters.noise_ps;
eps_min=parameters.eps_min;
beta=parameters.beta;  %ʽ3-65ʱ��ƽ������
psi=parameters.psi;
psi_min=parameters.psi_min;  %ʽ3-67��3-68��ֵ
psi_max=parameters.psi_max;
psi_pmin=parameters.psi_pmin;
psi_pmax=parameters.psi_pmax; 
psi_fr_old=parameters.psi_fr;
psi_peak=parameters.psi_peak;
qmax=parameters.qmax;
alpha=parameters.alpha; %ʽ3-63ƽ������
GH1=parameters.GH1;
gamma_old = parameters.gamma; %ǰһ֡���������
eps_old = parameters.eps;

Sf = zeros(len,1);
Sf(1)=ns_ps(1);   %�����������׵�Ƶ��ƽ��
Sf(end)=ns_ps(end);
for f=2:len-1
    Sf(f)=0.25*ns_ps(f-1)+0.5* ns_ps(f)+0.25* ns_ps(f+1);   %����������ʹ��3��hanning������ƽ��
end
% sig=sqrt(ns_ps);  %�����ź�Ƶ��
gamma=min(Sf./noise_ps,40);  % ��������ȣ���ֵ������������ֵΪ16dB����
% eps=alpha*X2_prev./noise_mu2 + (1-alpha)*max(gamma-1,0);   %��������ȣ�ʽ3-63   
eps=alpha*(GH1.^2).*gamma_old + (1-alpha)*max(gamma-1,0);
eps=max(eps_min,eps);  % ������������ȴ��� -25 dB����

A=eps./(1+eps);
v=A.*gamma;
ei_v=0.5*expint(v);
GH1=A.*exp(ei_v);    %����GH1����ʹ���������ƽ��������
% GH1=max(1,GH1); 

% ---���������������� P(Ho) ---------------

len2=len/2+1; 
psi=beta*psi+(1-beta)*eps_old(1:len2); %��ʽ3-65���������ƽ��       
psi_local=smoothing(psi,1);  %�Ӵ�ƽ������������ȵľֲ���ȫ��ƽ��ֵ
psi_global=smoothing(psi,5);  
C=log10(psi_max/psi_min); 

Plocal=zeros(len2,1);   % ���� P_local
imax=find(psi_local>=psi_max);
Plocal(imax)=1;
ibet=find(psi_local>psi_min & psi_local<psi_max);
Plocal(ibet)=log10(psi_local(ibet)/psi_min)/C;
    
    
Pglob=zeros(len2,1);   % ���� P_global
imax=find(psi_global>=psi_max);
Pglob(imax)=1;
ibet=find(psi_global>psi_min & psi_global<psi_max);
Pglob(ibet)=log10(psi_global(ibet)/psi_min)/C;
    
psi_fr=mean(psi);  % ���� Pframe
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

q2 = 1- Plocal.*Pglob*Pframe;  % ������������
% q2= min(qmax,q2); %��ֵ�������ʵ���qmax<1
q = [q2; flipud(q2(2:len2-1))]; %���Ʋ������γ�FFT�������� 

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

   
    


 



