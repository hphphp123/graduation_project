%一种双麦克语音增强基于迭代的噪声互功率谱估计。参考文献
%（An iterative noise cross-PSD estimation for two-microphone speech enhancement）

close all;clc;clear all;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\factory1.wav');
noise2=wavread('d:\noisex-92\factory2.wav');
[x1,noise11]=add_noisedata(s,noise1,fs,fs,-5); %第一路通道
[x2,noise12]=add_noisedata(s,noise2,fs,fs,-5); %第二路通道
%短时傅里叶变换
wlen=256;
SP=0.5;
shiftlen=0.5*wlen;
wnd=hamming(wlen);
x11=segment(x1,wlen,SP,wnd);    %分帧处理用到segment函数
x22=segment(x2,wlen,SP,wnd); 
framenum=size(x11,2);           %帧数
X1=fft(x11);
X2=fft(x22);

lamda1=0.9;
lamda2=0.6;
lamda3=0.8;
%第一帧
Px1x1(:,1)=X1(:,1).*conj(X1(:,1));
Px1x2(:,1)=X1(:,1).*conj(X2(:,1));
Px2x2(:,1)=X2(:,1).*conj(X2(:,1));
Pn1n2(:,1)=X1(:,1).*conj(X2(:,1));
H(:,1)=zeros(wlen,1);
for i=2:framenum
    Px1x1(:,i)=lamda2*Px1x1(:,i-1)+(1-lamda2)*(X1(:,i).*conj(X1(:,i)));
    Px1x2(:,i)=lamda2*Px1x2(:,i-1)+(1-lamda2)*(X1(:,i).*conj(X2(:,i)));
    Px2x2(:,i)=lamda2*Px2x2(:,i-1)+(1-lamda2)*(X2(:,i).*conj(X2(:,i)));
    Pn1n2(:,i)=lamda1*Pn1n2(:,i-1)+(1-lamda1)*(X1(:,i).*conj(X1(:,i))).*(1-H(:,i-1));
    cc=abs(X1(:,i).*conj(X2(:,i)));
    dd=abs(Pn1n2(:,i));
    Rpost(:,i)=max((cc./dd)-1,0);
    ccc=abs(X1(:,i-1).*conj(X2(:,i-1)));
    ddd=abs(Pn1n2(:,i-1));
    Rprio(:,i)=lamda3*(H(:,i-1).^2).*(ccc./ddd)+(1-lamda3)* Rpost(:,i);
    H(:,i)= Rprio(:,i)./(1+ Rprio(:,i)).*(abs(Px1x2(:,i))./(sqrt(Px1x1(:,i).*Px2x2(:,i))));
end
   
%得到相干滤波器之后，进行分别进行相干滤波，再重叠相加得到增强语音;
j=sqrt(-1);
Y=H.*abs(X1);
phase=angle(X1);
XXX=exp(j*phase);
Spec1=Y.*XXX;
y=zeros((framenum-1)*shiftlen+wlen,1);
for i=1:framenum
    start=(i-1)*shiftlen+1;
    spec1=Spec1(:,i);
    y(start:start+wlen-1)=y(start:start+wlen-1)...
        +real(ifft(spec1,wlen));
end
snr=SNR_singlech(s(1:N-1),y);fprintf(' snr=%5.1f\n',snr);
subplot(411),plot(s); title('clean');axis([1 N -1 1]);
subplot(412),plot(x1); title('x1'); axis([1 N -1 1]);
subplot(413),plot(x2); title('x2'); axis([1 N -1 1]);
subplot(414),plot(y); title('enhanced signal'); axis([1 N -1 1]);

sound (x1);
pause(1)
sound(y);
    
    
