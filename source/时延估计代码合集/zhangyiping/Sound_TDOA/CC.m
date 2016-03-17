clc
clear all
close all
%%
% *����������*

%--��Դ��ز���
fm=2000;    %Hz: ��Դ��Ƶ�ź����Ƶ�� ����0.5ms
Ts=0.2;     %s: ��Դ���� 0.2s
%--�������źŴ�����ز���
fs=10*fm;   %������ Ҳ����50us��һ����
ts=1/fs;    %�������
T=0.2;      %s: �Ѽ�����T�룬����һ��λ��
tMic=0:1/fs:T-1/fs;   %��������ʱ��
nMic=length(tMic);    %�������ݳ���
Rlen=nMic; %����صĳ���
%--�����������
v=340;      %m/s: ����
SNR=-10;    %dB
%--��Դ��MICλ��
Lco_S = [20,30];
Loc_A = [10,0];
Loc_B = [0,10];
Loc_C = [20,0];

%%
% *ģ����Դ�źż���MIC���յ����ź�*

%--������Դ
t = 0:ts:0.4; %������������ʱ��0.4�룬Ҳ������2������  
s = chirp(mod(t,0.2),0,0.2,fm,'linear');%Դ�ź�
%--������Դ�Ĳ���
figure();
plot((1:8001)*0.05,s);
xlabel('ʱ��/\itms');
ylabel('���');
title('�����ź�');

%--������Դ���MIC֮��ľ���
dis_SA = sqrt(sum((Lco_S-Loc_A).^2));
dis_SB = sqrt(sum((Lco_S-Loc_B).^2));
dis_SC = sqrt(sum((Lco_S-Loc_C).^2));
%--��Դ��MIC��ʱ���ӳ�
Lag_SA = dis_SA/v;
Lag_SB = dis_SB/v;
Lag_SC = dis_SC/v;
%--ת����������
diff_AB =round((Lag_SA-Lag_SB)/ts);
diff_AC =round((Lag_SA-Lag_SC)/ts);
%--MIC���յ�������
sigMicA=s(1:nMic);
sigMicB=s(1+diff_AB:nMic+diff_AB);
sigMicC=s(1+diff_AC:nMic+diff_AC);

sigMicA=awgn(sigMicA,SNR,'measured');
sigMicB=awgn(sigMicB,SNR,'measured');
sigMicC=awgn(sigMicC,SNR,'measured');
%--�ź�ʱ��ͼ
figure();subplot(3,1,1);
plot((1:4000)*0.05,sigMicA);
subplot(3,1,2);
plot((1:4000)*0.05,sigMicB);
subplot(3,1,3);
plot((1:4000)*0.05,sigMicC);

%%
% *��CC��Cross-Correlation���㷨����ʱ�Ӳ�*

%CC�㷨����ʱ��
rMicAB=xcorr(sigMicA,sigMicB,Rlen,'biased'); %��MIC A��B�źŻ����
rMicAC=xcorr(sigMicA,sigMicC,Rlen,'biased'); %��MIC A��C�źŻ����
            
[val,DelayDifferAB]=max(rMicAB);  %��������ֵ��λ���������ӳٲ
[val,DelayDifferAC]=max(rMicAC);

%�����ӳٲ����
delayDifferABRes=-(Rlen+1)+rMicAB(DelayDifferAB+(-3:3))*(DelayDifferAB+(-3:3))'/sum(rMicAB(DelayDifferAB+(-3:3)));
delayDifferACRes=-(Rlen+1)+rMicAC(DelayDifferAC+(-3:3))*(DelayDifferAC+(-3:3))'/sum(rMicAC(DelayDifferAC+(-3:3)));
distDiffAB=delayDifferABRes/fs*v;
distDiffAC=delayDifferACRes/fs*v;

%�������
funSrc=inline('[sqrt(sum((x-A).^2))-sqrt(sum((x-B).^2))-dAB;sqrt(sum((x-A).^2))-sqrt(sum((x-C).^2))-dAC]','x','A','B','C','dAB','dAC');
[Loc_reslut,err]=fsolve(funSrc,Loc_A,[],Loc_A,Loc_B,Loc_C,distDiffAB,distDiffAC);    %������ֵΪ����һ��MIC��λ��
%--��Դλ��ͼ
figure();
plot(Lco_S(1),Lco_S(2),'ro');
hold on;
plot(Loc_reslut(1),Loc_reslut(2),'o');
%--�������
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %�������
disp('Error=');
disp(Error);







