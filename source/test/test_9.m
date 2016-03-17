%论文9，本算法相对经典CF-GSC有两点改进，其一是将相干滤波器作为广义旁瓣抵消器的后置滤波器，
%克服了经典模式对相干滤波器的约束，提出了一种基于时间平滑参数的迭代噪声PSD估计算法，
%该方法利用时间递归平均技术，并结合迭代噪声PSD估计技术，获得的噪声PSD估计值更准确。
%基于改进的迭代再生PSD估计的相干滤波器设计.

close all;clc;clear all;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\factory1.wav');
noise2=wavread('d:\noisex-92\babble.wav');
[x1,noise11]=add_noisedata(s,noise1,fs,fs,0); %第一路通道
[x2,noise12]=add_noisedata(s,noise2,fs,fs,0); %第二路通道
%x1和x2都是时延对齐后的信号；
d=(x1+x2)/2;
n=x1-x2;

u=0.0001;
M=32;
w=zeros(M,1);
var=zeros(1,M);
yout=zeros(1,M);
%beta=0.0001;
for i=M:N
    input1=d(i);
    input2=n(i:-1:i-M+1);
    yout(i)=w(1:M)'*input2;
    ygsc(i)=input1-yout(i);
    var(i)=0.2*(n(i)^2)+0.8*var(i-1);
    %var(i)=n(i)^2+beta;
    U1=u/var(i);
    w=w+U1*ygsc(i)*input2;
end

wlen=256;
SP=0.5;
shiftlen=0.5*wlen;
wnd=hamming(wlen);
y1=segment(ygsc,wlen,SP,wnd);    %分帧处理用到segment函数
framenum=size(y1,2);           %帧数
Ygsc=fft(y1);

Hprio=zeros(wlen,framenum);
lamda=zeros(wlen,framenum);
gamma=zeros(wlen,framenum);
Rpost=zeros(wlen,framenum);
PY=zeros(wlen,framenum);
PY(:,1)=Ygsc(:,1).*conj(Ygsc(:,1));
PV=PY;
a=0.7;b=15;c=0.8;
%时间递归平均算法利用噪声一般对语音频谱具有不均匀影响这一特点，即有
%些频谱区域受噪声影响比其他一些区域所受影响更大。不同频谱的分量极有可能
%具有不同的实际 SNR，因此，在某一特定频带的实际 SNR 很低的时候，可以对噪
%声 PSD 按频带进行估计和更新。同样的，当某一频带存在语音的概率较低的时候，
%可以更新单个频带的噪声 PSD
%该算法中，
%噪声 PSD 估计是基于过去的噪声估计与当前带噪语音谱的加权平均，权重会根据
%每个频点的实际 SNR 来自适应改变。后验信噪比大，平滑参数大，避免过估计，反之，紧随着
%带噪语音的互功率谱变化。
for  k=2:framenum
    PY(:,k)=a*PY(:,k-1)+(1-a)*(Ygsc(:,k).*conj(Ygsc(:,k)));
    PV(:,k)= lamda(:,k-1).*PV(:,k-1)+(1-lamda(:,k-1)).*(Ygsc(:,k).*conj...
             (Ygsc(:,k))).*(1-Hprio(:,k-1));
    gamma(:,k)=(PY(:,k))./(PV(:,k-1)); 
    gamma(:,k) = min(gamma(:,k),100);
    m(:,k)=exp(-(b.*(gamma(:,k)-1.5)));
    lamda(:,k)=1./(1+m(:,k));
    ccc(:,k)=(PY(:,k))./(PV(:,k));
    ccc(:,k)=min(ccc(:,k),100);
    Rpost(:,k)=max((ccc(:,k)-1),0);
    Pprio(:,k)=c*(Hprio(:,k).^2).*ccc(:,k)+(1-c)* Rpost(:,k);
    Hprio(:,k)=  Pprio(:,k)./(1+  Pprio(:,k));
end

   
%得到相干滤波器之后，进行分别进行相干滤波，再重叠相加得到增强语音;
j=sqrt(-1);
Y=Hprio.*abs(Ygsc);
phase=angle(Ygsc);
XXX=exp(j*phase);
Spec1=Y.*XXX;
y=zeros((framenum-1)*shiftlen+wlen,1);
for i=1:framenum
    start=(i-1)*shiftlen+1;
    spec1=Spec1(:,i);
    y(start:start+wlen-1)=y(start:start+wlen-1)...
        +real(ifft(spec1,wlen));
end 

 snr1=SNR_singlech(s(1:N-1),y);fprintf(' snr1=%5.1f\n',snr1); 
 snr2=SNR_singlech(s,ygsc);fprintf(' snr2=%5.1f\n',snr2);
% 
subplot(411),plot(x1); title('x1');axis([1 N -1 1]);
subplot(412),plot(x2); title('x2'); axis([1 N -1 1]);
subplot(413),plot(ygsc);  title('ygsc');axis([1 N -1 1]);
subplot(414),plot(y); title('y'); axis([1 N -1 1]);
sound(x1);
pause(1)
sound(y);
%     
    
    
    

