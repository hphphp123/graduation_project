%基于麦克风小阵的多噪声环境语音增强算法，
%关键词：非平稳噪声：麦克风阵列：广义旁瓣抵消：相干滤波语音增强

close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\white.wav');
[x0,noise]=add_noisedata(s,noise1,fs,fs,5); %第一路添加白噪的带噪语音
xn=wavread('d:\语音文件\clean\sp21.wav'); %干扰噪声，假设目标信号已对齐，干扰信号有延迟
[x1,noise1]=add_noisedata(s,xn,fs,fs,5);%x1为第一路信号
k=15; %延迟点数
x1=x1';
% %r=zeros(1,2000);
r=wgn(1,N,-20);%产生-20dB高斯噪声
% % r=randn(1,2000); %noise
% % r=r/max(abs(r));
% noise2=wavread('d:\noisex-92\white.wav');
% [x1,noise]=add_noisedata(x1,noise1,fs,fs,5); %第二路添加白噪的带噪语音
xn=xn';
xn1=[r(1:k) xn(1:N-k)] ;          %delay signal
[x2,noise2]=add_noisedata(s,xn1,fs,fs,5);
n=1:N;
figure 
subplot(211)
plot(n,x1);
axis([1 N -1 1]);
title('模拟第一路麦克风信号');
subplot(212)
plot(n,x2);
axis([1 N -1 1]);
title('延迟后模拟第二路麦克风信号')
% y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
% x=find(y==max(y)) ;            %延迟
% delay_cc=x-length(x1)
% subplot(313)
% x2=x2';
% x3=[x2(delay:N) r(1:delay-1)];
% plot(n,x3);
% axis([1 N -1 1]);
% title('时延估计对齐后得到的第二路麦克风信号');
x1=x1';
% x2=x2';
s=s';
snr1=SNR_singlech(s(1:N),x1);fprintf(' 第一路麦克风信号snr1=%5.1f\n',snr1); 
snr2=SNR_singlech(s(1:N),x2);fprintf(' 第二路麦克风信号snr2=%5.1f\n',snr2);
sound(x1);
pause(1);
sound(x2);
% X1=fft(x1,length(x1));
% X2=fft(x2,length(x2));
% FXX=X1.*conj(X2);
% Amplitude=abs(X1).*abs(X2);
% GCC=FXX./Amplitude;
% R_XX=ifft(GCC);
% [val1,DelayDifferAB]=max(R_XX);
% delay_phat=length(x1)-DelayDifferAB+1
% X1=fft(x1,length(x1));
% X2=fft(x2,length(x2));
% FXX1=X1.*conj(X1);
% FXX2=X2.*conj(X2);
% Amplitude=sqrt(abs(FXX1.*FXX2));
% GCC=FXX./Amplitude;
% R_XX=ifft(GCC);
% [val1,DelayDifferAB]=max(R_XX);
% delay_scot=length(x1)-DelayDifferAB+1
% X1=fft(x1,length(x1));
% X2=fft(x2,length(x2));
% FXX1=X1.*conj(X1);
% FXX2=X2.*conj(X2);
% Amplitude=abs(FXX1);
% GCC=FXX./Amplitude;
% R_XX=ifft(GCC);
% [val1,DelayDifferAB]=max(R_XX);
% delay_roth=length(x1)-DelayDifferAB+1
% len=length(x1);
% x=[zeros(1000,1);x1];
% w=zeros(1000,1);%初始1000阶加权系数
% u=0.000005;
% y=zeros(len,1);
% for i=1:len-4000
%      y(i)=(x(i+999:-1:i))'*w;     
%      e=x2(i)-y(i); 
%      w=w+2*u*e*x(i+999:-1:i);
% end
% subplot(411),plot(x);grid on
% subplot(412),plot(x2);grid on
% subplot(413),plot(y);grid on
% subplot(414),plot(w);grid on
%wavwrite(y,'处理后信号');
% y1=xcorr(x,x2);
% %subplot(515),plot(y1);grid on
% %delay_lms=find(y1==max(y1))-len
% delay_lms=find(w==max(w))
%%%%%%%%%%%%%%%%%模拟通道信号完成%%%%%%%%%%%%%%
wlen=256;
SP=0.5;
shiftlen=0.5*wlen;
wnd=hamming(wlen);
y1=segment(x1,wlen,SP,wnd);    %分帧处理用到segment函数
y2=segment(x2,wlen,SP,wnd);
framenum=size(y1,2);           %帧数
Y1=fft(y1);
Y2=fft(y2);
%
% 
for i=1:framenum
    for j=1:wlen
        Pxy(j,i)=Y1(j,i).*conj(Y2(j,i));
        Pxx(j,i)=Y1(j,i).*conj(Y1(j,i));
        Pyy(j,i)=Y2(j,i).*conj(Y2(j,i));
       Txy(j,i)=real(Pxy(j,i)./sqrt(Pxx(j,i).*Pyy(j,i)));
    end
end

T1=0.35;T2=0.8;
for i=1:framenum
    for j=1:wlen
          if Txy(j,i)<T1
             Txy(j,i)=0.001;
         elseif Txy(j,i)>T2
             Txy(j,i)=0.99;
          end
    end
end

%得到相干滤波器Txy之后，进行分别进行相干滤波，再重叠相加得到新的x1和x2;
j=sqrt(-1);
Y11=Txy.*abs(Y1);
Y22=Txy.*abs(Y2);
phase=angle(Y1);
XXX=exp(j*phase);
Spec1=Y11.*XXX;
Spec2=Y22.*XXX;
y11=zeros((framenum-1)*shiftlen+wlen,1);
y22=zeros((framenum-1)*shiftlen+wlen,1);
for i=1:framenum
    start=(i-1)*shiftlen+1;
    spec1=Spec1(:,i);
    spec2=Spec2(:,i);
    y11(start:start+wlen-1)=y11(start:start+wlen-1)...
        +real(ifft(spec1,wlen));
    y22(start:start+wlen-1)=y22(start:start+wlen-1)...
        +real(ifft(spec2,wlen));
end

%广义旁瓣抵消语音增强部分算法，包括固定波束形成和噪声估计模块。

d=(y11+y22)/2;     %固定波束形成
%d1=(x1+x2)/2;
%在固定波束形成后加谱减法
%d1=simplesubspec(signal,wlen,inc,NIS,a,b);
n=y11-y22;         %噪声信号
%n1=x1-x2;
%自适应噪声抵消模块，使用归一化最小均方算法迭代进行噪声估计；

u=0.001;
M=32;
w=zeros(M,1);
var=zeros(1,M);
%yout=zeros(1,N-1);
for i=M:N-1
    input1=d(i);
    input2=n(i:-1:i-M+1);
    yout(i)=w(1:M)'*input2;
    e(i)=input1-yout(i);
    var(i)=0.2*n(i)^2+0.8*var(i-1);
    U1=u/var(i);
    w=w+U1.*e(i).*input2;
end
%u1=0.0001;var1=zeros(1,M);
% for i=M:N
%     input1=d1(i);
%     input2=n1(i:-1:i-M+1);
%     yout(i)=w(1:M)*input2;
%     e1(i)=input1-yout(i);
%     var1=0.2*(input2)'.^2+0.8*var1;
%     U1=u./var1;
%     w=w+U1.*(e1(i).*input2');
% end
% 

 snr3=SNR_singlech(s(1:N-1),d);fprintf(' 相干滤波后的信号snr3=%5.1f\n',snr3); 
 snr4=SNR_singlech(s(1:N-1),e);fprintf(' 经过GSC后信号snr4=%5.1f\n',snr4);
%  snr3=SNR_singlech(s,d1);fprintf(' snr3=%5.1f\n',snr3); 
%  snr4=SNR_singlech(s,e1);fprintf(' snr4=%5.1f\n',snr4);
n=1:length(s);
figure 
subplot(411)
plot(n+1,x1);
axis([1 length(s)-1 -1.5 1.5]);
title('第一路语音信号');
subplot(412)
plot(n+1,x2);
axis([1 length(s)-1 -1.5 1.5]);
title('第二路语音信号');
n=1:length(s)-1;
subplot(413)
plot(n,d);
axis([1 length(s)-1 -1.5 1.5]);
title('固定波束形成后信号');
subplot(414)
plot(n,e);
axis([1 length(s)-1 -1.5 1.5]);
title('经过相干滤波和tf-GSC后的信号');


 sound(d);
 pause(1)
 sound(e);


