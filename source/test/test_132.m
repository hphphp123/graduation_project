%������˷�С��Ķ���������������ǿ�㷨��
%�ؼ��ʣ���ƽ����������˷����У������԰����������˲�������ǿ
close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\�����ļ�\clean\sp01.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
xn=wavread('d:\�����ļ�\clean\sp21.wav'); %��������������Ŀ���ź��Ѷ��룬�����ź����ӳ�
[x1,noise1]=add_noisedata(s,xn,fs,fs,5);%x1Ϊ��һ·�ź�
k=15; %�ӳٵ���
x1=x1';
r=wgn(1,N,-20);%����-20dB��˹����
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
title('ԭʼ�����ź�');
subplot(312)
plot(n,x1);
axis([1 N -1 1]);
title('ģ���һ·��˷��ź�');
subplot(313)
plot(n,x2);
axis([1 N -1 1]);
title('�ӳٺ�ģ��ڶ�·��˷��ź�')
plot(n,x3);
axis([1 N -1 1]);
title('�ӳٺ�ģ�����·��˷��ź�')
%%%%%%%%%%%%%%%%%ģ��ͨ���ź����%%%%%%%%%%%%%%
wlen=200;
inc=80;
win=hanning(wlen);
N=length(x1);time=(0:N-1)/fs;
y=enframe(x,win,inc)';

fn=size(y1,2);           %֡��
frameTime=(((1:fn)-1)*80+wlen/2)/fs; %ÿ֡��Ӧʱ��
W2=wlen/2+1;n2=1:W2;freq=(n2-1)*fs/wlen;
Y1=fft(x1);
Y2=fft(x2);
Y3=fft(x3);
figure
set(gcf,'Position',[20 100 600 500]);
axes('Position',[0.1 0.1 0.85 0.5]);
imagesc(frameTime,freq,abs(Y1(n2)));
axis xy;ylabel('Ƶ��/Hz');xlabel('ʱ��/s');
title('����ͼ');
m=64;
LightYellow=[0.6 0.6 0.6];
MidRed=[0 0 0];
Black=[0.5 0.7 1];
Colors=[LightYellow MidRed Black];
colormap(SpecColorMap(m,Colors));
