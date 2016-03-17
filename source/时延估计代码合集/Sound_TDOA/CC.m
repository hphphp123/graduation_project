clc
clear all
close all
%%
% *各参数设置*

%--声源相关参数
fm=2000;    %Hz: 信源调频信号最高频率 周期0.5ms
Ts=0.2;     %s: 信源周期 0.2s
%--采样和信号处理相关参数
fs=10*fm;   %采样率 也就是50us采一次样
ts=1/fs;    %采样间隔
T=0.2;      %s: 搜集数据T秒，计算一次位置
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
Rlen=nMic; %做相关的长度
%--物理参数设置
v=340;      %m/s: 音速
SNR=-10;    %dB
%--声源和MIC位置
Lco_S = [20,30];
Loc_A = [10,0];
Loc_B = [0,10];
Loc_C = [20,0];

%%
% *模拟声源信号及各MIC接收到的信号*

%--产生声源
t = 0:ts:0.4; %假设声波持续时间0.4秒，也就是有2个周期  
s = chirp(mod(t,0.2),0,0.2,fm,'linear');%源信号
%--画出声源的波形
figure();
plot((1:8001)*0.05,s);
xlabel('时间/\itms');
ylabel('振幅');
title('声音信号');

%--计算信源与个MIC之间的距离
dis_SA = sqrt(sum((Lco_S-Loc_A).^2));
dis_SB = sqrt(sum((Lco_S-Loc_B).^2));
dis_SC = sqrt(sum((Lco_S-Loc_C).^2));
%--信源到MIC的时间延迟
Lag_SA = dis_SA/v;
Lag_SB = dis_SB/v;
Lag_SC = dis_SC/v;
%--转化成相差点数
diff_AB =round((Lag_SA-Lag_SB)/ts);
diff_AC =round((Lag_SA-Lag_SC)/ts);
%--MIC接收到的数据
sigMicA=s(1:nMic);
sigMicB=s(1+diff_AB:nMic+diff_AB);
sigMicC=s(1+diff_AC:nMic+diff_AC);

sigMicA=awgn(sigMicA,SNR,'measured');
sigMicB=awgn(sigMicB,SNR,'measured');
sigMicC=awgn(sigMicC,SNR,'measured');
%--信号时域图
figure();subplot(3,1,1);
plot((1:4000)*0.05,sigMicA);
subplot(3,1,2);
plot((1:4000)*0.05,sigMicB);
subplot(3,1,3);
plot((1:4000)*0.05,sigMicC);

%%
% *用CC（Cross-Correlation）算法估计时延差*

%CC算法求延时差
rMicAB=xcorr(sigMicA,sigMicB,Rlen,'biased'); %求MIC A、B信号互相关
rMicAC=xcorr(sigMicA,sigMicC,Rlen,'biased'); %求MIC A、C信号互相关
            
[val,DelayDifferAB]=max(rMicAB);  %互相关最大值的位置体现了延迟差。
[val,DelayDifferAC]=max(rMicAC);

%最终延迟差估计
delayDifferABRes=-(Rlen+1)+rMicAB(DelayDifferAB+(-3:3))*(DelayDifferAB+(-3:3))'/sum(rMicAB(DelayDifferAB+(-3:3)));
delayDifferACRes=-(Rlen+1)+rMicAC(DelayDifferAC+(-3:3))*(DelayDifferAC+(-3:3))'/sum(rMicAC(DelayDifferAC+(-3:3)));
distDiffAB=delayDifferABRes/fs*v;
distDiffAC=delayDifferACRes/fs*v;

%方程求解
funSrc=inline('[sqrt(sum((x-A).^2))-sqrt(sum((x-B).^2))-dAB;sqrt(sum((x-A).^2))-sqrt(sum((x-C).^2))-dAC]','x','A','B','C','dAB','dAC');
[Loc_reslut,err]=fsolve(funSrc,Loc_A,[],Loc_A,Loc_B,Loc_C,distDiffAB,distDiffAC);    %迭代初值为其中一个MIC的位置
%--声源位置图
figure();
plot(Lco_S(1),Lco_S(2),'ro');
hold on;
plot(Loc_reslut(1),Loc_reslut(2),'o');
%--计算误差
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %估计误差
disp('Error=');
disp(Error);







