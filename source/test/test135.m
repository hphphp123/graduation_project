%固定波束形成路加vad检测和维纳滤波，消除非相干噪声

close all;clc;clear all;
[s,fs,bits]=wavread('bluesky3.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\white.wav');
[x0,noise]=add_noisedata(s,noise1,fs,fs,100); %第一路添加白噪的带噪语音


%模拟另一个方向上的干扰信号，每个有一个相等的k个点的延迟，初始信干比为5dB，主瓣方向语音信号已对齐
k=5; %延迟点数
r=wgn(1,N,-20); %产生-20dB高斯噪声
xn=wavread('doct3.wav'); %干扰噪声，假设目标信号已对齐，干扰信号有延迟
xn=xn-mean(xn);
xn=xn/max(abs(xn));
xn=xn';
[x1,noise1]=add_noisedata(x0,xn,fs,fs,0); %第一路添加干扰的带噪语音
xn1=[r(1:k) xn(1:N-k)];
xn1=xn1';
[x2,noise2]=add_noisedata(x0,xn1,fs,fs,0); %第二路添加延迟干扰的带噪语音
xn2=[r(1:2*k) xn(1:N-2*k)];
xn2=xn2';
[x3,noise3]=add_noisedata(x0,xn2,fs,fs,0); %第二路添加延迟干扰的带噪语音
xn3=[r(1:3*k) xn(1:N-3*k)];
xn3=xn3';
[x4,noise4]=add_noisedata(x0,xn3,fs,fs,0); %第二路添加延迟干扰的带噪语音
% sound(x0,fs);
% sound(xn,fs);
% sound(x1,fs);
% sound(x2,fs);
% sound(x3,fs);
% sound(x4,fs);

%画出四路波形图
% figure(1);
% time=1:N;
% subplot 211 ;
% plot(time, x0(time));
% subplot 212;
% plot(time, xn(time));
% figure(2);
% time=1:N;
% subplot 411;
% plot(time,x1(time));
% subplot 412;
% plot(time,x2(time));
% subplot 413;
% plot(time,x3(time));
% subplot 414;
% plot(time,x4(time));
M=4;
Len = length(x1);
mic1 = x1(1:Len);
mic2 = x2(1:Len);
mic3 = x3(1:Len);
mic4 = x4(1:Len);

x = [mic1';mic2';mic3';mic4'];
FBFout=sum(x)/M;

B = [1 -1 0 0 
     0 1 -1 0 
     0 0 1 -1 ];
Bout = B*x;
x1 = Bout(1,:)';
x2 = Bout(2,:)';
x3 = Bout(3,:)';
%x4 = Bout(4,:)';

y1 = zeros(1,Len);
y2 = y1;
y3 = y1;
y4 = y1;
MCout = y1;
N = 64;
h1 = zeros(1,N);
h2 = h1;
h3 = h1; 
h4 = h1;
GSCout = zeros(1,Len);
u = 0.0002;
weight = [];
for i=1:fix(Len/N)-1
    X1 = fft(x1((i-1)*N+1:(i+1)*N));
    H1 = fft([h1,zeros(1,N)]);
    O11 = real(ifft(X1'.*H1));
    y1(i*N+1:(i+1)*N) = O11(N+1:N*2);
    
    X2 = fft(x2((i-1)*N+1:(i+1)*N));
    H2 = fft([h2,zeros(1,N)]);
    O12 = real(ifft(X2'.*H2));
    y2(i*N+1:(i+1)*N) = O12(N+1:N*2);
    
    X3 = fft(x3((i-1)*N+1:(i+1)*N));
    H3 = fft([h3,zeros(1,N)]);
    O13 = real(ifft(X3'.*H3));    
    y3(i*N+1:(i+1)*N) = O13(N+1:N*2);
    
%     X4 = fft(x4((i-1)*N+1:(i+1)*N));
%     H4 = fft([h4,zeros(1,N)]);
%     O14 = real(ifft(X4'.*H4));
%     y4(i*N+1:(i+1)*N) = O14(N+1:N*2);
    
    X = [X1';X2';X3'];%;X4'
    MCin = sum(X); 
    MCout(i*N+1:(i+1)*N) = sum([y1(i*N+1:(i+1)*N);y2(i*N+1:(i+1)*N);y3(i*N+1:(i+1)*N)]);%;y4(i*N+1:(i+1)*N)
    
    GSCout(i*N+1:(i+1)*N) = FBFout(i*N+1:(i+1)*N)-(MCout(i*N+1:(i+1)*N));
    E = fft([zeros(1,N),GSCout(i*N+1:(i+1)*N)]);
    O2 = real(ifft(E.*conj(MCin)));
    V = O2(1:N);
    
    h1 = h1+2*u*V;
    h2 = h1+2*u*V;
    h3 = h3+2*u*V;
    %h4 = h4+2*u*V;
    
    weight = [weight h1(1)];
end
subplot(211);plot(mic1);sound(mic1)
subplot(212);plot(GSCout);sound(GSCout)
snr2=SNR_singlech(s,mic1);fprintf('snr2=%5.1f\n',snr2);
snr1=SNR_singlech(s,GSCout);fprintf('snr1=%5.1f\n',snr1);
%wavwrite(GSCout,'Fgsc.wav');


    

    
    
    
    
