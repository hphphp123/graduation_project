%基于麦克风小阵的多噪声环境语音增强算法，
%关键词：非平稳噪声：麦克风阵列：广义旁瓣抵消：相干滤波语音增强
close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
xn=wavread('d:\语音文件\clean\sp21.wav'); %干扰噪声，假设目标信号已对齐，干扰信号有延迟
[x1,noise1]=add_noisedata(s,xn,fs,fs,5);%x1为第一路信号
k=15; %延迟点数
x1=x1';
r=wgn(1,N,-20);%产生-20dB高斯噪声
xn=xn';
xn1=[r(1:k) xn(1:N-k)] ;          %delay signal
[x2,noise2]=add_noisedata(s,xn1,fs,fs,5);
xn=xn';
xn2=[r(1:k) xn1(1:N-k)] ;          %delay signal
[x3,noise3]=add_noisedata(s,xn2,fs,fs,5);
n=1:N;
figure 
subplot(311)
plot(n,s);
title('原始语音信号');
subplot(312)
plot(n,x1);
axis([1 N -1 1]);
title('模拟第一路麦克风信号');
subplot(313)
plot(n,x2);
axis([1 N -1 1]);
title('延迟后模拟第二路麦克风信号')
plot(n,x3);
axis([1 N -1 1]);
title('延迟后模拟第三路麦克风信号')
%%%%%%%%%%%%%%%%%模拟通道信号完成%%%%%%%%%%%%%%
wlen=200;
inc=80;
win=hanning(wlen);
N=length(x1);time=(0:N-1)/fs;
y=enframe(x,win,inc)';

fn=size(y1,2);           %帧数
frameTime=(((1:fn)-1)*80+wlen/2)/fs; %每帧对应时间
W2=wlen/2+1;n2=1:W2;freq=(n2-1)*fs/wlen;
Y1=fft(x1);
Y2=fft(x2);
Y3=fft(x3);
figure
set(gcf,'Position',[20 100 600 500]);
axes('Position',[0.1 0.1 0.85 0.5]);
imagesc(frameTime,freq,abs(Y1(n2)));
axis xy;ylabel('频率/Hz');xlabel('时间/s');
title('语谱图');
m=64;
LightYellow=[0.6 0.6 0.6];
MidRed=[0 0 0];
Black=[0.5 0.7 1];
Colors=[LightYellow MidRed Black];
colormap(SpecColorMap(m,Colors));
