%����9�����㷨��Ծ���CF-GSC������Ľ�����һ�ǽ�����˲�����Ϊ�����԰�������ĺ����˲�����
%�˷��˾���ģʽ������˲�����Լ���������һ�ֻ���ʱ��ƽ�������ĵ�������PSD�����㷨��
%�÷�������ʱ��ݹ�ƽ������������ϵ�������PSD���Ƽ�������õ�����PSD����ֵ��׼ȷ��
%���ڸĽ��ĵ�������PSD���Ƶ�����˲������.

close all;clc;clear all;
[s,fs,bits]=wavread('d:\�����ļ�\clean\sp01.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
noise1=wavread('d:\noisex-92\factory1.wav');
noise2=wavread('d:\noisex-92\babble.wav');
[x1,noise11]=add_noisedata(s,noise1,fs,fs,0); %��һ·ͨ��
[x2,noise12]=add_noisedata(s,noise2,fs,fs,0); %�ڶ�·ͨ��
%x1��x2����ʱ�Ӷ������źţ�
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
y1=segment(ygsc,wlen,SP,wnd);    %��֡�����õ�segment����
framenum=size(y1,2);           %֡��
Ygsc=fft(y1);

Hprio=zeros(wlen,framenum);
lamda=zeros(wlen,framenum);
gamma=zeros(wlen,framenum);
Rpost=zeros(wlen,framenum);
PY=zeros(wlen,framenum);
PY(:,1)=Ygsc(:,1).*conj(Ygsc(:,1));
PV=PY;
a=0.7;b=15;c=0.8;
%ʱ��ݹ�ƽ���㷨��������һ�������Ƶ�׾��в�����Ӱ����һ�ص㣬����
%ЩƵ������������Ӱ�������һЩ��������Ӱ����󡣲�ͬƵ�׵ķ������п���
%���в�ͬ��ʵ�� SNR����ˣ���ĳһ�ض�Ƶ����ʵ�� SNR �ܵ͵�ʱ�򣬿��Զ���
%�� PSD ��Ƶ�����й��ƺ͸��¡�ͬ���ģ���ĳһƵ�����������ĸ��ʽϵ͵�ʱ��
%���Ը��µ���Ƶ�������� PSD
%���㷨�У�
%���� PSD �����ǻ��ڹ�ȥ�����������뵱ǰ���������׵ļ�Ȩƽ����Ȩ�ػ����
%ÿ��Ƶ���ʵ�� SNR ������Ӧ�ı䡣��������ȴ�ƽ�������󣬱�������ƣ���֮��������
%���������Ļ������ױ仯��
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

   
%�õ�����˲���֮�󣬽��зֱ��������˲������ص���ӵõ���ǿ����;
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
    
    
    

