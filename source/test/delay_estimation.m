close all;clc;clear all;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\white.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,5); %第一路添加白噪的带噪语音
k=1000;
x1=x1';
%r=zeros(1,2000);
r=wgn(1,N,-20);%产生-20dB高斯噪声
% r=randn(1,2000); %noise
% r=r/max(abs(r));
x2=[r(1:k) x1(1:N-k)] ;          %delay signal
% noise2=wavread('d:\noisex-92\babble.wav');
% [x2,noise]=add_noisedata(x2,noise2,fs,fs,5); %第二路添加白噪的带噪语音
n=1:N;
figure 
subplot(311)
plot(n,x1);
axis([1 N -1 1]);
title('模拟第一路麦克风信号');
subplot(312)
plot(n,x2);
axis([1 N -1 1]);
title('延迟后模拟第二路麦克风信号')
y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %延迟
delay=x-length(x1)
subplot(313)
%x2=x2';
x3=[x2(delay:N) r(1:delay-1)];
plot(n,x3);
axis([1 N -1 1]);
title('时延估计对齐后得到的第二路麦克风信号');
x1=x1';
x2=x2';
s=s';