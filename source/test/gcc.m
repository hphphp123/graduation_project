close all;clc;clear all;
[s,fs,bits]=wavread('d:\�����ļ�\clean\sp01.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
noise1=wavread('d:\noisex-92\white.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,5); %��һ·���Ӱ���Ĵ�������
k=1000;
x1=x1';

%r=zeros(1,2000);
r=wgn(1,N,-20);%����-20dB��˹����
% r=randn(1,2000); %noise
% r=r/max(abs(r));
x2=[r(1:k) x1(1:N-k)] ;          %delay signal
% noise2=wavread('d:\noisex-92\babble.wav');
% [x2,noise]=add_noisedata(x2,noise2,fs,fs,5); %�ڶ�·���Ӱ���Ĵ�������
% wlen=256;
% SP=0.5;
% shiftlen=0.5*wlen;
% wnd=hamming(wlen);
% y1=segment(x1,wlen,SP,wnd);    %��֡�����õ�segment����
% y2=segment(x2,wlen,SP,wnd);
% framenum=size(y1,2);           %֡��

% for i=1:framenum
%     for j=1:wlen
%         Pxy(j,i)=Y1(j,i).*conj(Y2(j,i));
%         Pxx(j,i)=Y1(j,i).*conj(Y1(j,i));
%         Pyy(j,i)=Y2(j,i).*conj(Y2(j,i));
%     end
% end

n=1:N;
figure 
subplot(311)
plot(n,x1);
axis([1 N -1 1]);
title('ģ���һ·��˷��ź�');
subplot(312)
plot(n,x2);
axis([1 N -1 1]);
title('�ӳٺ�ģ��ڶ�·��˷��ź�')
y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %�ӳ�
delay1=x-length(x1)
subplot(313)
%x2=x2';
%x3=[x2(delay:N) r(1:delay-1)];
plot(y);
title('ʱ�ӹ��ƶ����õ��ĵڶ�·��˷��ź�');
