%������˷�С��Ķ���������������ǿ�㷨��
%�ؼ��ʣ���ƽ����������˷����У������԰����������˲�������ǿ

close all;clc;clear all;
[s,fs,bits]=wavread('d:\�����ļ�\clean\sp01.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
noise1=wavread('d:\noisex-92\white.wav');
noise2=wavread('d:\noisex-92\babble.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,-5); %��һ·��Ӱ���Ĵ�������
[x2,noise]=add_noisedata(s,noise2,fs,fs,100); %�ڶ�·��Ӱ���Ĵ�������
%x1��x2����ʱ�Ӷ������źţ�
%
%������ͨ������ɺ���
%���źŽ��з�֡�����������Ϊ8khz,֡��256��֡��128��ÿ֡�Ӻ����������ж�ʱ����Ҷ�任
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

    

    
    
    
    
