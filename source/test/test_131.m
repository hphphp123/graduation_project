%������˷�С��Ķ���������������ǿ�㷨��
%�ؼ��ʣ���ƽ����������˷����У������԰����������˲�������ǿ

close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\�����ļ�\clean\sp01.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
noise1=wavread('d:\noisex-92\white.wav');
[x0,noise]=add_noisedata(s,noise1,fs,fs,5); %��һ·��Ӱ���Ĵ�������
xn=wavread('d:\�����ļ�\clean\sp21.wav'); %��������������Ŀ���ź��Ѷ��룬�����ź����ӳ�
[x1,noise1]=add_noisedata(s,xn,fs,fs,5);%x1Ϊ��һ·�ź�
k=15; %�ӳٵ���
x1=x1';
% %r=zeros(1,2000);
r=wgn(1,N,-20);%����-20dB��˹����
% % r=randn(1,2000); %noise
% % r=r/max(abs(r));
% noise2=wavread('d:\noisex-92\white.wav');
% [x1,noise]=add_noisedata(x1,noise1,fs,fs,5); %�ڶ�·��Ӱ���Ĵ�������
xn=xn';
xn1=[r(1:k) xn(1:N-k)] ;          %delay signal
[x2,noise2]=add_noisedata(s,xn1,fs,fs,5);
n=1:N;
figure 
subplot(211)
plot(n,x1);
axis([1 N -1 1]);
title('ģ���һ·��˷��ź�');
subplot(212)
plot(n,x2);
axis([1 N -1 1]);
title('�ӳٺ�ģ��ڶ�·��˷��ź�')
% y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
% x=find(y==max(y)) ;            %�ӳ�
% delay_cc=x-length(x1)
% subplot(313)
% x2=x2';
% x3=[x2(delay:N) r(1:delay-1)];
% plot(n,x3);
% axis([1 N -1 1]);
% title('ʱ�ӹ��ƶ����õ��ĵڶ�·��˷��ź�');
x1=x1';
% x2=x2';
s=s';
snr1=SNR_singlech(s(1:N),x1);fprintf(' ��һ·��˷��ź�snr1=%5.1f\n',snr1); 
snr2=SNR_singlech(s(1:N),x2);fprintf(' �ڶ�·��˷��ź�snr2=%5.1f\n',snr2);
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
% w=zeros(1000,1);%��ʼ1000�׼�Ȩϵ��
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
%wavwrite(y,'������ź�');
% y1=xcorr(x,x2);
% %subplot(515),plot(y1);grid on
% %delay_lms=find(y1==max(y1))-len
% delay_lms=find(w==max(w))
%%%%%%%%%%%%%%%%%ģ��ͨ���ź����%%%%%%%%%%%%%%
wlen=256;
SP=0.5;
shiftlen=0.5*wlen;
wnd=hamming(wlen);
y1=segment(x1,wlen,SP,wnd);    %��֡�����õ�segment����
y2=segment(x2,wlen,SP,wnd);
framenum=size(y1,2);           %֡��
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

%�õ�����˲���Txy֮�󣬽��зֱ��������˲������ص���ӵõ��µ�x1��x2;
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

%�����԰����������ǿ�����㷨�������̶������γɺ���������ģ�顣

d=(y11+y22)/2;     %�̶������γ�
%d1=(x1+x2)/2;
%�ڹ̶������γɺ���׼���
%d1=simplesubspec(signal,wlen,inc,NIS,a,b);
n=y11-y22;         %�����ź�
%n1=x1-x2;
%����Ӧ��������ģ�飬ʹ�ù�һ����С�����㷨���������������ƣ�

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

 snr3=SNR_singlech(s(1:N-1),d);fprintf(' ����˲�����ź�snr3=%5.1f\n',snr3); 
 snr4=SNR_singlech(s(1:N-1),e);fprintf(' ����GSC���ź�snr4=%5.1f\n',snr4);
%  snr3=SNR_singlech(s,d1);fprintf(' snr3=%5.1f\n',snr3); 
%  snr4=SNR_singlech(s,e1);fprintf(' snr4=%5.1f\n',snr4);
n=1:length(s);
figure 
subplot(411)
plot(n+1,x1);
axis([1 length(s)-1 -1.5 1.5]);
title('��һ·�����ź�');
subplot(412)
plot(n+1,x2);
axis([1 length(s)-1 -1.5 1.5]);
title('�ڶ�·�����ź�');
n=1:length(s)-1;
subplot(413)
plot(n,d);
axis([1 length(s)-1 -1.5 1.5]);
title('�̶������γɺ��ź�');
subplot(414)
plot(n,e);
axis([1 length(s)-1 -1.5 1.5]);
title('��������˲���tf-GSC����ź�');


 sound(d);
 pause(1)
 sound(e);


