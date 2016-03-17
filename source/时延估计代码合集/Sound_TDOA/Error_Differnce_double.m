clc
clear all
close all
err_data=[];
%%
% *����������*
for  x=0:1:100 
    
%--��Դ��ز���
fmin=500;
fmax=2000;    %Hz: ��ԴΪһƵ�ʽ���������źţ����Ƶ��fmin�����Ƶ��fmax
S_last=0.2;   %s ����Դ����ʱ��
%--�������źŴ�����ز���
fs=3e6;   %������ 
ts=1/fs;    %�������

T=0.3;      %s: �Ѽ�����T�룬����һ��λ��
tMic=0:1/fs:T-1/fs;   %��������ʱ��
nMic=length(tMic);    %�������ݳ���
%--�����������
v=34000;      %cm/s: ����
%--��Դ��MICλ��   ��λCM
d=20;    %��˷���
Lco_S = [50,x];
Loc_A = [0,0];
Loc_B = [d,0];
Loc_C = [2*d,0];

Loc_E = [100,100];
Loc_F = [100-d,100];
Loc_G = [100-2*d,100];
%%
% *ģ����Դ�źż���MIC���յ����ź�*

%--������Դ
t = 0:ts:S_last; %������������ʱ��0.2�� 
s = chirp(t,fmin,S_last,fmax,'linear');%Դ�ź�, Ƶ�����Ե����������ź�
nsource=length(s);

%--������Դ������MIC֮����ӳٵ���
diff_SA = sqrt(sum((Lco_S-Loc_A).^2))/v/ts;
diff_SB = sqrt(sum((Lco_S-Loc_B).^2))/v/ts;
diff_SC = sqrt(sum((Lco_S-Loc_C).^2))/v/ts;

diff_SE = sqrt(sum((Lco_S-Loc_E).^2))/v/ts;
diff_SF = sqrt(sum((Lco_S-Loc_F).^2))/v/ts;
diff_SG = sqrt(sum((Lco_S-Loc_G).^2))/v/ts;

%--MIC���յ�������
sigMicA=zeros(1,nMic);
sigMicB=zeros(1,nMic);
sigMicC=zeros(1,nMic);

sigMicA(1+diff_SA:nsource+diff_SA)=s;
sigMicB(1+diff_SB:nsource+diff_SB)=s;
sigMicC(1+diff_SC:nsource+diff_SC)=s;

sigMicE=zeros(1,nMic);
sigMicF=zeros(1,nMic);
sigMicG=zeros(1,nMic);

sigMicE(1+diff_SE:nsource+diff_SE)=s;
sigMicF(1+diff_SF:nsource+diff_SF)=s;
sigMicG(1+diff_SG:nsource+diff_SG)=s;

%%
% *��GCC��Generalized Cross-Correlation�� Phase Transform �㷨����ʱ�Ӳ�*   

% ��һ�鴫����

% Phase Transform�㷨����ʱ��
FA=fft(sigMicA,length(sigMicA));  %����Ҷ�任
FB=fft(sigMicB,length(sigMicB));
FAB=FA.*conj(FB);               %�󻥹�����
Amplitude=abs(FA).*abs(FB);     %�����
GCC=FAB./Amplitude;             %ȥ��������Ϣ
R_ab=ifft(GCC);                 
[val1,DelayDifferAB]=max(R_ab);  %��������ֵ��λ���������ӳٲ

FC=fft(sigMicC,length(sigMicC));
FAC=FA.*conj(FC);               %�󻥹�����
Amplitude=abs(FA).*abs(FC);     %�����
GCC=FAC./Amplitude;             %ȥ��������Ϣ
R_ac=ifft(GCC);                 
[val2,DelayDifferAC]=max(R_ac);  %��������ֵ��λ���������ӳٲ

if(DelayDifferAB<3000)
  distDiffAB=(DelayDifferAB-1)/fs*v;
else distDiffAB=(DelayDifferAB-nMic-1)/fs*v;
end

if(DelayDifferAC<3000)
  distDiffAC=(DelayDifferAC-1)/fs*v;
else distDiffAC=(DelayDifferAC-nMic-1)/fs*v;
end

flag=0;
for a=10:20:100
    for b=10:20:100
%�������
funSrc=inline('[sqrt(sum((x-A).^2))-sqrt(sum((x-B).^2))-dAB;sqrt(sum((x-A).^2))-sqrt(sum((x-C).^2))-dAC]','x','A','B','C','dAB','dAC');
[Loc_reslut1,err]=fsolve(funSrc,[a,b],[100,100],Loc_A,Loc_B,Loc_C,distDiffAB,distDiffAC);    %������ֵΪ����һ��MIC��λ��
Error1=sqrt(sum((Loc_reslut1-Lco_S).^2)); %�������
 if Error1<1;
    flag=1;
    break;    
 end
    end
    if flag ==1;
      break;
    end
end

%%
% *��GGG��Generalized Gross-Gorrelation�� Phase Transform �㷨����ʱ�Ӳ�*

% �ڶ��鴫����

% Phase Transform�㷨����ʱ��
FE=fft(sigMicE,length(sigMicE));  %����Ҷ�任
FF=fft(sigMicF,length(sigMicF));
FEF=FE.*conj(FF);               %�󻥹�����
amplitude=abs(FE).*abs(FF);     %�����
GCC=FEF./amplitude;             %ȥ��������Ϣ
R_ab=ifft(GCC);                 
[val1,DelayDifferEF]=max(R_ab);  %��������ֵ��λ���������ӳٲ

FG=fft(sigMicG,length(sigMicG));
FEG=FE.*conj(FG);               %�󻥹�����
amplitude=abs(FE).*abs(FG);     %�����
GCC=FEG./amplitude;             %ȥ��������Ϣ
R_ac=ifft(GCC);                 
[val2,DelayDifferEG]=max(R_ac);  %��������ֵ��λ���������ӳٲ

if(DelayDifferEF<3000)
  distDiffEF=(DelayDifferEF-1)/fs*v;
else distDiffEF=(DelayDifferEF-nMic-1)/fs*v;
end

if(DelayDifferEG<3000)
  distDiffEG=(DelayDifferEG-1)/fs*v;
else distDiffEG=(DelayDifferEG-nMic-1)/fs*v;
end

flag=0;
for a=10:20:100
    for b=10:20:100
%�������
funSrc=inline('[sqrt(sum((x-E).^2))-sqrt(sum((x-F).^2))-dEF;sqrt(sum((x-E).^2))-sqrt(sum((x-G).^2))-dEG]','x','E','F','G','dEF','dEG');
[Loc_reslut2,err]=fsolve(funSrc,[a,b],[100,100],Loc_E,Loc_F,Loc_G,distDiffEF,distDiffEG);    %������ֵΪ����һ��MIG��λ��
Error2=sqrt(sum((Loc_reslut2-Lco_S).^2)); %�������
 if Error2<1;
    flag=1;
    break;    
 end
    end
    if flag ==1;
      break;
    end
end
%%
% *ȡ����ȷ�ľ���*

if ( Error1<1 &&  Error2<1)
   Loc_reslut=(Loc_reslut1+Loc_reslut2)/2;
else if (Error2<1)
    Loc_reslut=Loc_reslut2;
     else Loc_reslut=Loc_reslut1;
    end
end 

%--�������
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %�������
err_data=[err_data Error];

end
       figure();
       plot(0:100,err_data,'-*');
       ylabel('���')
       xlabel('������50�µ�������')
       grid on;
       box on;
       
          
set(gcf,'Units','centimeters','Position',[10 10 6.67 5]);%
set(gcf,'Units','centimeters','Position',[10 10 9 6.75]);%
%set(gcf,'Units','centimeters','Position',[10 10 13.5 9]);%

set(get(gca,'XLabel'),'FontSize',8);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',8);
set(get(gca,'TITLE'),'FontSize',8);
set(gca,'fontsize',8);
set(gca,'linewidth',0.5); %�����ߴ�0.5��
%set(gca,'box','off');%Controls the box around the plotting area
set(get(gca,'Children'),'linewidth',1.2);%����ͼ���߿�1.5��















