clc;clear;
%设置输入设备参数
% ai = analoginput('mcc');
% addchannel(ai,0:2);
% set(ai,'samplerate',fs);
% T=0.5;      %s: 搜集数据T秒，计算一次位置
% buffer = fs * T; %单次采样时间，对应采样点数
% set(ai,'SamplesPerTrigger',buffer);
% start(ai);
fs=100000;   %采样率
nMic=50000;
v=340000;
d=60;
% data = getdata(ai);
% Phase Transform算法求延时差
load original.mat
sigMicA_in = data(1:nMic,1);
sigMicB_in = data(1:nMic,2);
sigMicC_in = data(1:nMic,3);
x1=sigMicA_in;
x2=sigMicB_in;
x3=sigMicC_in;
Loc_A = [-d,0];
Loc_B = [0,0];
Loc_C = [d,0];
%%  广义互相关cc
N=length(x1);
n=1:N;
y=xcorr(x1,x3);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %延迟
%delay=x-length(x1);
delay=x;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_cc1=delay/fs*v


N=length(x1);
n=1:N;
y=xcorr(x3,x2);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %延迟
%delay=x-length(x1);
delay=x;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_cc2=delay/fs*v


A=delay_cc2;
B=delay_cc1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;
%     if()

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta_cc = acos(cosA)*180/pi;



%% 广义互相关 phat
X1=fft(x1,length(x1));
X3=fft(x3,length(x3));
FXX=X1.*conj(X3);
Amplitude=abs(FXX);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay=length(x1)-DelayDifferAB+1;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_phat1=delay/fs*v

X2=fft(x2,length(x2));
X3=fft(x3,length(x3));
FXX=X2.*conj(X3);
Amplitude=abs(FXX);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay=length(x1)-DelayDifferAB+1;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_phat2=delay/fs*v

A=delay_phat2;
B=delay_phat1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta_phat = acos(cosA)*180/pi;

%% 广义互相关 scot
X1=fft(x1,length(x1));
X3=fft(x3,length(x3));
FXX1=X1.*conj(X1);
FXX3=X3.*conj(X3);
FXX=X1.*conj(X3);
Amplitude=sqrt(abs(FXX1.*FXX3))+0.1;
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay=length(x1)-DelayDifferAB+1;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_scot1=delay/fs*v

X2=fft(x2,length(x2));
X3=fft(x3,length(x3));
FXX2=X2.*conj(X2);
FXX3=X3.*conj(X3);
FXX=X2.*conj(X3);
Amplitude=sqrt(abs(FXX2.*FXX3));
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay=length(x1)-DelayDifferAB+1;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_scot2=delay/fs*v


A=delay_scot2;
B=delay_scot1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta_scot = acos(cosA)*180/pi;
%% 广义互相关 roth
X1=fft(x1,length(x1));
X3=fft(x3,length(x3));
FXX1=X1.*conj(X1);
FXX3=X3.*conj(X3);
Amplitude=abs(FXX3);
FXX=X1.*conj(X3);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay=length(x1)-DelayDifferAB+1;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_roth1=delay/fs*v

X2=fft(x2,length(x2));
X3=fft(x3,length(x3));
FXX2=X2.*conj(X2);
FXX3=X3.*conj(X3);
FXX=X2.*conj(X3);
Amplitude=abs(FXX3);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay=length(x1)-DelayDifferAB+1;
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_roth2=delay/fs*v

A=delay_roth2;
B=delay_roth1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta_roth = acos(cosA)*180/pi;
%% lms
len=length(x1);
x=[zeros(1000,1);x1];
w=zeros(1000,1);%初始1000阶加权系数，滤波器长度设为1000
u=0.000005;
y=zeros(len,1);
for i=1:len-4000
     y(i)=(x(i+999:-1:i))'*w;     
     e=x3(i)-y(i); 
     w=w+2*u*e*x(i+999:-1:i);
end

delay=find(w==max(w));
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_lms1=delay/fs*v

len=length(x2);
x=[zeros(1000,1);x2];
w=zeros(1000,1);%初始1000阶加权系数，滤波器长度设为1000
u=0.000005;
y=zeros(len,1);
for i=1:len-4000
     y(i)=(x(i+999:-1:i))'*w;     
     e=x3(i)-y(i); 
     w=w+2*u*e*x(i+999:-1:i);
end

delay=find(w==max(w));
if(delay<5000)
  delay=(delay);
else delay=(delay-nMic);
end
delay_lms2=delay/fs*v

A=delay_lms2;
B=delay_lms1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta_lms = acos(cosA)*180/pi;