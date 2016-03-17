clc
clear all
close all
%%
% *各参数设置*

%--声源相关参数
fmin=500;
fmax=2000;    %Hz: 信源为一频率渐变的余弦信号，最低频率fmin，最高频率fmax
S_last=0.2;   %s ：声源持续时间
%--采样和信号处理相关参数
fs=3e6;   %采样率 
ts=1/fs;    %采样间隔

T=0.3;      %s: 搜集数据T秒，计算一次位置
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
%--物理参数设置
v=34000;      %cm/s: 音速
%--声源和MIG位置   单位GM
d=10;    %麦克风间距
Lco_S = [80,20];
Loc_E = [100,100];
Loc_F = [100-d,100];
Loc_G = [100-2*d,100];
%%
% *模拟声源信号及各MIG接收到的信号*

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
%--计算信源传到各MIG之间的延迟点数
diff_SE = sqrt(sum((Lco_S-Loc_E).^2))/v/ts;
diff_SF = sqrt(sum((Lco_S-Loc_F).^2))/v/ts;
diff_SG = sqrt(sum((Lco_S-Loc_G).^2))/v/ts;

% %--计算信源与个MIG之间的距离
% dis_SE = sqrt(sum((Lco_S-Loc_E).^2));
% dis_SF = sqrt(sum((Lco_S-Loc_F).^2));
% dis_SG = sqrt(sum((Lco_S-Loc_G).^2));
% %--信源到MIG的时间延迟
% Lag_SE = dis_SE/v;
% Lag_SF = dis_SF/v;
% Lag_SG = dis_SG/v;
% %--转化成相差点数
% diff_EF =round((Lag_SE-Lag_SF)/ts);
% diff_EG =round((Lag_SE-Lag_SG)/ts);

%--MIG接收到的数据
sigMicE=zeros(1,nMic);
sigMicF=zeros(1,nMic);
sigMicG=zeros(1,nMic);

sigMicE(1+diff_SE:nsource+diff_SE)=s;
sigMicF(1+diff_SF:nsource+diff_SF)=s;
sigMicG(1+diff_SG:nsource+diff_SG)=s;

%--信号时域图
figure();subplot(3,1,1);
plot(sigMicE);
subplot(3,1,2);
plot(sigMicF);
subplot(3,1,3);
plot(sigMicG);
%%
% *用GGG（Generalized Gross-Gorrelation） Phase Transform 算法估计时延差*

% Phase Transform算法求延时差
FE=fft(sigMicE,length(sigMicE));  %傅里叶变换
FF=fft(sigMicF,length(sigMicF));
FEF=FE.*conj(FF);               %求互功率谱
amplitude=abs(FE).*abs(FF);     %求幅度
GCC=FEF./amplitude;             %去除幅度信息
R_ab=ifft(GCC);                 
[val1,DelayDifferEF]=max(R_ab);  %互相关最大值的位置体现了延迟差。

FG=fft(sigMicG,length(sigMicG));
FEG=FE.*conj(FG);               %求互功率谱
amplitude=abs(FE).*abs(FG);     %求幅度
GCC=FEG./amplitude;             %去除幅度信息
R_ac=ifft(GCC);                 
[val2,DelayDifferEG]=max(R_ac);  %互相关最大值的位置体现了延迟差。

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
%方程求解
funSrc=inline('[sqrt(sum((x-E).^2))-sqrt(sum((x-F).^2))-dEF;sqrt(sum((x-E).^2))-sqrt(sum((x-G).^2))-dEG]','x','E','F','G','dEF','dEG');
[Loc_reslut,err]=fsolve(funSrc,[a,b],[100,100],Loc_E,Loc_F,Loc_G,distDiffEF,distDiffEG);    %迭代初值为其中一个MIG的位置
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %估计误差
 if Error<1;
    flag=1;
    break;    
 end
    end
    if flag ==1;
      break;
    end
end
%--声源位置图
figure();
plot(Lco_S(1),Lco_S(2),'ro');
hold on;
plot(Loc_reslut(1),Loc_reslut(2),'o');
%--计算误差
Frror=sqrt(sum((Loc_reslut-Lco_S).^2)); %估计误差
disp('Error=');
disp(Error);













