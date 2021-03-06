clc
clear all
close all
%%
% *各参数设置*

%--声源相关参数
fmin=500;
fmax=2000;    %Hz: 信源为一频率渐变的余弦信号，最低频率fmin，最高频率fmax
S_last=0.1;   %s ：声源持续时间
%--采样和信号处理相关参数
fs=3e4;   %采样率 
ts=1/fs;    %采样间隔

T=0.12;      %s: 搜集数据T秒，计算一次位置
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
%--物理参数设置
v=34000;      %cm/s: 音速
%--声源和MIC位置   单位CM
d=5;    %麦克风间距
Lco_S = [10,10];
Loc_A = [50-d,0];
Loc_B = [50,0];
Loc_C = [50+d,0];
%%
% *模拟声源信号及各MIC接收到的信号*

%--产生声源
t = 0:ts:S_last; %假设声波持续时间0.2秒 
s = chirp(t,fmin,S_last,fmax,'linear');%源信号, 频率线性递增的余弦信号
nsource=length(s);
%--画出声源的波形
figure();
plot(s);
xlabel('时间/\itus');
ylabel('振幅');
title('声音信号');

%--计算信源传到各MIC之间的延迟点数
diff_SA = sqrt(sum((Lco_S-Loc_A).^2))/v/ts;
diff_SB = sqrt(sum((Lco_S-Loc_B).^2))/v/ts;
diff_SC = sqrt(sum((Lco_S-Loc_C).^2))/v/ts;

%--计算信源与个MIC之间的距离
dis_SA = sqrt(sum((Lco_S-Loc_A).^2));
dis_SB = sqrt(sum((Lco_S-Loc_B).^2));
dis_SC = sqrt(sum((Lco_S-Loc_C).^2));
%--信源到MIC的时间延迟
dis_AC=dis_SA-dis_SC;
dis_BC=dis_SB-dis_SC;

Lag_SA = dis_SA/v;
Lag_SB = dis_SB/v;
Lag_SC = dis_SC/v;
%--转化成相差点数
diff_AC =round((Lag_SA-Lag_SC)/ts);
diff_BC =round((Lag_SB-Lag_SC)/ts);

%--MIC接收到的数据
sigMicA=zeros(1,nMic);
sigMicB=zeros(1,nMic);
sigMicC=zeros(1,nMic);

sigMicA(1+diff_SA:nsource+diff_SA)=s;
sigMicB(1+diff_SB:nsource+diff_SB)=s;
sigMicC(1+diff_SC:nsource+diff_SC)=s;

%--信号时域图
figure();subplot(3,1,1);
plot(sigMicA);
subplot(3,1,2);
plot(sigMicB);
subplot(3,1,3);
plot(sigMicC);
%%
% *用GCC（Generalized Cross-Correlation） Phase Transform 算法估计时延差*

% Phase Transform算法求延时差
FB=fft(sigMicB,length(sigMicB));  %傅里叶变换
FC=fft(sigMicC,length(sigMicC));
FBC=FB.*conj(FC);               %求互功率谱
Amplitude=abs(FB).*abs(FC);     %求幅度
GCC=FBC./Amplitude;             %去除幅度信息
R_BC=ifft(GCC);                 
[val1,DelayDifferBC]=max(R_BC);  %互相关最大值的位置体现了延迟差。

FA=fft(sigMicA,length(sigMicA));  %傅里叶变换
FC=fft(sigMicC,length(sigMicC));
FAC=FA.*conj(FC);               %求互功率谱
Amplitude=abs(FA).*abs(FC);     %求幅度
GCC=FAC./Amplitude;             %去除幅度信息
R_AC=ifft(GCC);                 
[val1,DelayDifferAC]=max(R_AC);  %互相关最大值的位置体现了延迟差。

if(DelayDifferBC<50)
  distDiffBC=(DelayDifferBC-1)/fs*v;
else distDiffBC=(DelayDifferBC-nMic-1)/fs*v;
end

if(DelayDifferAC<50)
  distDiffAC=(DelayDifferAC-1)/fs*v;
else distDiffAC=(DelayDifferAC-nMic-1)/fs*v;
end

A=distDiffBC;
B=distDiffAC;

M=2*A^2-B^3+6*A*B^2-10*A^2*B+4*d^2*A-2*d^2*B;
N=-8*A^2-2*B^2+8*A*B;
r3=M/N;


r2=M1/N1;

cos=(r3^2+d^2-r2^2)/(2*d*r3);
sin=sqrt(1-cos^2);

a=(50-d)+r3*cos;
b=r3*sin;
Loc_reslut=[a,b];
%--声源位置图
figure();
plot(Lco_S(1),Lco_S(2),'ro');
hold on;
plot(a,b,'o');
%--计算误差
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %估计误差
disp('Error=');
disp(Error);













