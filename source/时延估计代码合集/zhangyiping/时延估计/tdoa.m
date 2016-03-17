%������˷�С��Ķ���������������ǿ�㷨��
%�ؼ��ʣ���ƽ����������˷����У������԰����������˲�������ǿ

close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\�����ļ�\clean\sp01.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
noise1=wavread('d:\noisex-92\white.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,5); %��һ·��Ӱ���Ĵ�������
k=5;
x1=x1';
%r=zeros(1,2000);
r=wgn(1,N,-20);%����-20dB��˹����
% r=randn(1,2000); %noise
% r=r/max(abs(r));
x2=[r(1:k) x1(1:N-k)] ;          %delay signal
noise2=wavread('d:\noisex-92\babble.wav');
[x2,noise]=add_noisedata(x2,noise2,fs,fs,5); %�ڶ�·��Ӱ���Ĵ�������
%%  ���廥���cc
n=1:N;
y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %�ӳ�
delay_cc=x-length(x1)
%% ���廥��� phat
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
%% ���廥��� scot
X1=fft(x1,length(x1));
X2=fft(x2,length(x2));
FXX1=X1.*conj(X1);
FXX2=X2.*conj(X2);
Amplitude=sqrt(abs(FXX1.*FXX2));
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_scot=length(x1)-DelayDifferAB+1
%% ���廥��� roth
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
w=zeros(1000,1);%��ʼ1000�׼�Ȩϵ�����˲���������Ϊ1000
u=0.000005;
y=zeros(len,1);
for i=1:len-4000
     y(i)=(x(i+999:-1:i))'*w;     
     e=x2(i)-y(i); 
     w=w+2*u*e*x(i+999:-1:i);
end

delay_lms=find(w==max(w))