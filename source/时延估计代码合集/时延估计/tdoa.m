%基于麦克风小阵的多噪声环境语音增强算法，
%关键词：非平稳噪声：麦克风阵列：广义旁瓣抵消：相干滤波语音增强

close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\white.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,5); %第一路添加白噪的带噪语音
k=5;
x1=x1';
%r=zeros(1,2000);
r=wgn(1,N,-20);%产生-20dB高斯噪声
% r=randn(1,2000); %noise
% r=r/max(abs(r));
x2=[r(1:k) x1(1:N-k)] ;          %delay signal
noise2=wavread('d:\noisex-92\babble.wav');
[x2,noise]=add_noisedata(x2,noise2,fs,fs,5); %第二路添加白噪的带噪语音
%%  广义互相关cc
n=1:N;
y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %延迟
delay_cc=x-length(x1)
%% 广义互相关 phat
x1=x1';
s=s';
X1=fft(x1,length(x1));
X2=fft(x2,length(x2));
FXX=X1.*conj(X2);
Amplitude=abs(X1).*abs(X2);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_phat=length(x1)-DelayDifferAB+1
%% 广义互相关 scot
X1=fft(x1,length(x1));
X2=fft(x2,length(x2));
FXX1=X1.*conj(X1);
FXX2=X2.*conj(X2);
Amplitude=sqrt(abs(FXX1.*FXX2));
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_scot=length(x1)-DelayDifferAB+1
%% 广义互相关 roth
X1=fft(x1,length(x1));
X2=fft(x2,length(x2));
FXX1=X1.*conj(X1);
FXX2=X2.*conj(X2);
Amplitude=abs(FXX1);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_roth=length(x1)-DelayDifferAB+1
%% lms
len=length(x1);
x=[zeros(1000,1);x1];
w=zeros(1000,1);%初始1000阶加权系数，滤波器长度设为1000
u=0.000005;
y=zeros(len,1);
for i=1:len-4000
     y(i)=(x(i+999:-1:i))'*w;     
     e=x2(i)-y(i); 
     w=w+2*u*e*x(i+999:-1:i);
end

delay_lms=find(w==max(w))