%基于麦克风小阵的多噪声环境语音增强算法，
%关键词：非平稳噪声：麦克风阵列：广义旁瓣抵消：相干滤波语音增强

close all;clc;clear all;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\white.wav');
noise2=wavread('d:\noisex-92\babble.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,-5); %第一路添加白噪的带噪语音
[x2,noise]=add_noisedata(s,noise2,fs,fs,100); %第二路添加白噪的带噪语音
%x1和x2都是时延对齐后的信号；
%
%计算两通道的相干函数
%对信号进行分帧，这里采样率为8khz,帧长256，帧移128，每帧加汉明窗，进行短时傅里叶变换
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

 snr1=SNR_singlech(s(1:N-1),d);fprintf(' snr1=%5.1f\n',snr1); 
 snr2=SNR_singlech(s(1:N-1),e);fprintf(' snr2=%5.1f\n',snr2);
%  snr3=SNR_singlech(s,d1);fprintf(' snr3=%5.1f\n',snr3); 
%  snr4=SNR_singlech(s,e1);fprintf(' snr4=%5.1f\n',snr4);

 sound(y11);
 pause(1)
 sound(e);
%  pause(1);
%  sound(e1);
% %  
figure 
subplot(511),plot(s);title('clean signal');  axis([1 N -1 1]);
subplot(512),plot(x1);title('first signal');   axis([1 N -1 1]);
subplot(513),plot(x2);title('second signal');  axis([1 N -1 1]);
subplot(514),plot(e);title('enhanced signal'); axis([1 N -1 1]);
subplot(515),plot(d);title('enhanced signal1'); axis([1 N -1 1]);
% figure
% subplot(411),plot(x1);title('x1');
% subplot(412),plot(x2);title('x2');    axis([1 N -1 1]);
% subplot(413),plot(y11);title('y11');   axis([1 N -1 1]);
% subplot(414),plot(y22);title('y22');  axis([1 N -1 1]);
%  figure
% subplot(411),plot(y11);title('y11');  axis([1 N -1 1]);
% subplot(412),plot(y22);title('y22');  axis([1 N -1 1]);
% subplot(413),plot(d);title('d');  axis([1 N -1 1]);
% subplot(414),plot(e);title('e');  axis([1 N -1 1]);
%     
%     

    

    
    
    
    
