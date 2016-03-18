%�̶������γ�·��vad����ά���˲����������������

close all;clc;clear all;
[s,fs,bits]=wavread('bluesky3.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
% noise1=wavread('d:\noisex-92\white.wav');
snr=0;
% [x0,noise]=add_noisedata(s,noise1,fs,fs,snr); %��һ·��Ӱ���Ĵ�������
x0=s;
%ģ����һ�������ϵĸ����źţ�ÿ����һ����ȵ�k������ӳ٣���ʼ�Ÿɱ�Ϊ5dB�����귽�������ź��Ѷ���
k=5; %�ӳٵ���
r=wgn(1,N,-20); %����-20dB��˹����
xn=wavread('doct3.wav'); %��������������Ŀ���ź��Ѷ��룬�����ź����ӳ�
xn=xn-mean(xn);
xn=xn/max(abs(xn));
xn=xn';
[x1,noise1]=add_noisedata(x0,xn,fs,fs,snr); %��һ·��Ӹ��ŵĴ�������
xn1=[r(1:k) xn(1:N-k)];
xn1=xn1';
[x2,noise2]=add_noisedata(x0,xn1,fs,fs,snr); %�ڶ�·����ӳٸ��ŵĴ�������
xn2=[r(1:2*k) xn(1:N-2*k)];
xn2=xn2';
[x3,noise3]=add_noisedata(x0,xn2,fs,fs,snr); %�ڶ�·����ӳٸ��ŵĴ�������
xn3=[r(1:3*k) xn(1:N-3*k)];
xn3=xn3';
[x4,noise4]=add_noisedata(x0,xn3,fs,fs,snr); %�ڶ�·����ӳٸ��ŵĴ�������

k1=3;
r=wgn(1,N,-20); %����-20dB��˹����
xn=wavread('d:\noisex-92\f16.wav'); %��������������Ŀ���ź��Ѷ��룬�����ź����ӳ�
xn=xn-mean(xn);
xn=xn/max(abs(xn));
xn=xn';
[x1,noise1]=add_noisedata(x1,xn,fs,fs,snr); %��һ·��Ӹ��ŵĴ�������
xn1=[r(1:k1) xn(1:N-k1)];
xn1=xn1';
[x2,noise2]=add_noisedata(x2,xn1,fs,fs,snr); %�ڶ�·����ӳٸ��ŵĴ�������
xn2=[r(1:2*k1) xn(1:N-2*k1)];
xn2=xn2';
[x3,noise3]=add_noisedata(x3,xn2,fs,fs,snr); %�ڶ�·����ӳٸ��ŵĴ�������
xn3=[r(1:3*k1) xn(1:N-3*k1)];
xn3=xn3';
[x4,noise4]=add_noisedata(x4,xn3,fs,fs,snr); %�ڶ�·����ӳٸ��ŵĴ�������
% sound(x0,fs);
% wavwrite(x0,fs,'x0');
% sound(xn,fs);
% sound(x1,fs);
% sound(x2,fs);
% sound(x3,fs);
% sound(x4,fs);

%������·����ͼ
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

% ������ͨ������ɺ���
% ���źŽ��з�֡�����������Ϊ8khz,֡��256��֡��128��ÿ֡�Ӻ����������ж�ʱ����Ҷ�任
wlen=256;
SP=0.5;
shiftlen=0.5*wlen;
wnd=hamming(wlen);
y1=segment(x1,wlen,SP,wnd);    %��֡�����õ�segment����
y2=segment(x3,wlen,SP,wnd);
framenum=size(y1,2);           %֡��
Y1=fft(y1);
Y2=fft(y2);
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
y3=segment(x2,wlen,SP,wnd);    %��֡�����õ�segment����
y4=segment(x4,wlen,SP,wnd);
framenum=size(y1,2);           %֡��
Y3=fft(y3);
Y4=fft(y4);
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


j=sqrt(-1);
Y33=Txy.*abs(Y3);
Y44=Txy.*abs(Y4);
phase=angle(Y3);
XXX=exp(j*phase);
Spec1=Y33.*XXX;
Spec2=Y44.*XXX;
y33=zeros((framenum-1)*shiftlen+wlen,1);
y44=zeros((framenum-1)*shiftlen+wlen,1);
for i=1:framenum
    start=(i-1)*shiftlen+1;
    spec1=Spec1(:,i);
    spec2=Spec2(:,i);
    y33(start:start+wlen-1)=y33(start:start+wlen-1)...
        +real(ifft(spec1,wlen));
    y44(start:start+wlen-1)=y44(start:start+wlen-1)...
        +real(ifft(spec2,wlen));
end
%�����԰����������ǿ�����㷨�������̶������γɺ���������ģ�顣
d=(y11+y22+y33+y44)/4;
% d=(x1+x2+x3+x4)/4;     %�̶������γ�
% sound(d);
% figure(3);
% time=1:N;
% plot(time,d);
%d1=(x1+x2)/2;
% snr0=SNR_singlech(s,x0);            % ����ά���˲���������
% snr1=SNR_singlech(s,x1)
% snr2=SNR_singlech(s,x2) 
% snr3=SNR_singlech(s,x3) 
% snr4=SNR_singlech(s,x4) 
% snr5=SNR_singlech(s,d) 


% 
%�ڹ̶������γ�һ·���ά���˲���
% IS=.15; %�޻���ʱ��
% output=WienerScalart96m_2(d,fs,IS,0.12);
% d=output;

%����Ӧ����ģ�飬�˲�������Ϊ512
B=[1,-1,0,0;0,1,-1,0;0,0,1,-1]
u=0.001;
sysorder=32;
N=4; 
M=length(d);
w1=zeros(sysorder,1);
w2=zeros(sysorder,1);
w3=zeros(sysorder,1);
y=zeros(1,sysorder);
W=[w1';w2';w3'];
for n=sysorder+1:M
    X1=x1(n-sysorder+1:1:n);
    X2=x2(n-sysorder+1:1:n);
    X3=x3(n-sysorder+1:1:n);
    X4=x4(n-sysorder+1:1:n);
    X=[X1';X2';X3';X4'];
    Z=B*X;
    Y=W.*Z;
    y(n)=sum(Y(:));
    e(n)=d(n)-y(n);
    W=W+u*e(n)*Z;
end
% figure(1)
% subplot 211;
% plot(s);title('Ŀ�������ź�');
% sound(s);pause(2);
% subplot 212;
% plot(xn(1:length(d)));title('�����ź�');
% sound(xn(1:length(d)));pause(2);
% figure(2)
% subplot 211;
% plot(x1);
% sound(x1);title('��һ·�ɼ����������ź�');pause(2);axis([0 2e4 min(x1) max(x1) ]);
% subplot 212;
% plot(d);title('�̶������γɺ��ź�');
% sound(d);pause(2);
% figure(3)
% subplot 211;
% plot(y);title('���Ƴ����ĸ����ź�');axis([0 2e4 -0.5 0.5]);
% sound(y);pause(2);
% subplot 212;
% plot(e);title('����Ӧ�����γɺ��ź�');axis([0 2e4 -1 1]);
% sound(e);pause(2);
% m=4;
% Wc=[0.25;0.25;0.25;0.25];
% wop=Wc;
% drawpp(m,wop);
% wop=B'*W(:,sysorder);
% drawpp(m,wop);
% wop=Wc-B'*W(:,sysorder);
% drawpp(m,wop);
% % theta=[0:180];
% % [theta,f]=meshgrid(theta,f);
% % figure(4)
% % mesh(theta,f,20*log10(abs(W)));
% name1=strcat('y1','.pcm');
% fid=fopen(name1,'w');
% fwrite(fid,y(1,:)*32767/3,'int16');
% name2=strcat('out_p','.pcm');
% fid=fopen(name2,'w');
% fwrite(fid,e(1,:)*32767/3,'int16');
% name3=strcat('xn','.pcm');
% fid=fopen(name3,'w');
% fwrite(fid,xn(1,:)*32767/3,'int16');
% fclose all;
snr0=SNR_singlech(s(1:M),e)           % ����ά���˲���������
wavwrite(e,fs,'f16-0-0-cf');
snr1=SNR_singlech(s(1:M),d)
snr2=SNR_singlech(s,x1)