clc
clear all
close all
%%
% *����������*

%--��Դ��ز���
fmin=500;
fmax=2000;    %Hz: ��ԴΪһƵ�ʽ���������źţ����Ƶ��fmin�����Ƶ��fmax
S_last=0.1;   %s ����Դ����ʱ��
%--�������źŴ�����ز���
fs=3e6;   %������ 
ts=1/fs;    %�������

T=0.12;      %s: �Ѽ�����T�룬����һ��λ��
tMic=0:1/fs:T-1/fs;   %��������ʱ��
nMic=length(tMic);    %�������ݳ���
%--�����������
v=34000;      %cm/s: ����
%--��Դ��MICλ��   ��λCM
d=30;    %��˷���
Lco_S = [10,10];
Loc_A = [50-d,0];
Loc_B = [50,0];
Loc_C = [50+d,0];
%%
% *ģ����Դ�źż���MIC���յ����ź�*

%--������Դ
t = 0:ts:S_last; %������������ʱ��0.2�� 
s = chirp(t,fmin,S_last,fmax,'linear');%Դ�ź�, Ƶ�����Ե����������ź�
nsource=length(s);
%--������Դ�Ĳ���
figure();
plot(s);
xlabel('ʱ��/\itus');
ylabel('���');
title('�����ź�');

%--������Դ������MIC֮����ӳٵ���
diff_SA = sqrt(sum((Lco_S-Loc_A).^2))/v/ts;
diff_SB = sqrt(sum((Lco_S-Loc_B).^2))/v/ts;
diff_SC = sqrt(sum((Lco_S-Loc_C).^2))/v/ts;

%--������Դ���MIC֮��ľ���
dis_SA = sqrt(sum((Lco_S-Loc_A).^2));
dis_SB = sqrt(sum((Lco_S-Loc_B).^2));
dis_SC = sqrt(sum((Lco_S-Loc_C).^2));
%--��Դ��MIC��ʱ���ӳ�
dis_AC=dis_SA-dis_SC;
dis_BC=dis_SB-dis_SC;

Lag_SA = dis_SA/v;
Lag_SB = dis_SB/v;
Lag_SC = dis_SC/v;
%--ת����������
diff_AC =round((Lag_SA-Lag_SC)/ts);
diff_BC =round((Lag_SB-Lag_SC)/ts);

%--MIC���յ�������
sigMicA=zeros(1,nMic);
sigMicB=zeros(1,nMic);
sigMicC=zeros(1,nMic);

sigMicA(1+diff_SA:nsource+diff_SA)=s;
sigMicB(1+diff_SB:nsource+diff_SB)=s;
sigMicC(1+diff_SC:nsource+diff_SC)=s;

%--�ź�ʱ��ͼ
figure();subplot(3,1,1);
plot(sigMicA);
subplot(3,1,2);
plot(sigMicB);
subplot(3,1,3);
plot(sigMicC);
%%
% *��GCC��Generalized Cross-Correlation�� Phase Transform �㷨����ʱ�Ӳ�*

% Phase Transform�㷨����ʱ��
FB=fft(sigMicB,length(sigMicB));  %����Ҷ�任
FC=fft(sigMicC,length(sigMicC));
FBC=FB.*conj(FC);               %�󻥹�����
Amplitude=abs(FB).*abs(FC);     %�����
GCC=FBC./Amplitude;             %ȥ��������Ϣ
R_BC=ifft(GCC);                 
[val1,DelayDifferBC]=max(R_BC);  %��������ֵ��λ���������ӳٲ

FA=fft(sigMicA,length(sigMicA));  %����Ҷ�任
FC=fft(sigMicC,length(sigMicC));
FAC=FA.*conj(FC);               %�󻥹�����
Amplitude=abs(FA).*abs(FC);     %�����
GCC=FAC./Amplitude;             %ȥ��������Ϣ
R_AC=ifft(GCC);                 
[val1,DelayDifferAC]=max(R_AC);  %��������ֵ��λ���������ӳٲ

if(DelayDifferBC<5000)
  distDiffBC=(DelayDifferBC-1)/fs*v;
else distDiffBC=(DelayDifferBC-nMic-1)/fs*v;
end

if(DelayDifferAC<5000)
  distDiffAC=(DelayDifferAC-1)/fs*v;
else distDiffAC=(DelayDifferAC-nMic-1)/fs*v;
end

A=distDiffBC;
B=distDiffAC;

M=4*A^3-B^3+6*A*B^2-10*A^2*B+4*d^2*A-2*d^2*B;
N=-8*A^2-2*B^2+8*A*B;
r3=M/N;

M1=4*A^3-B^3+4*A*B^2-6*A^2*B+2*d^2*B-4*d^2*A;
N1=8*A^2+2*B^2-8*A*B;
r2=M1/N1;

cos=(r3^2+d^2-r2^2)/(2*d*r3);
sin=sqrt(1-cos^2);

a=(50-d)+r3*cos;
b=r3*sin;
Loc_reslut=[a,b];
%--��Դλ��ͼ
figure();
plot(Lco_S(1),Lco_S(2),'ro');
hold on;
plot(a,b,'o');
%--�������
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %�������
disp('Error=');
disp(Error);













