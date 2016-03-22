%�̶������γ�·��vad����ά���˲����������������

close all;clc;clear all;
[s,fs,bits]=wavread('tone4.wav');                 %�������ź�
s=s-mean(s);                   %����ֱ������
s=s/max(abs(s));               %��ֵ��һ��
N=length(s);                   %�����źų���
% noise1=wavread('d:\noisex-92\white.wav');
snr=-5;
% [x0,noise]=add_noisedata(s,noise1,fs,fs,snr); %��һ·��Ӱ���Ĵ�������
x0=s;

%ģ����һ�������ϵĸ����źţ�ÿ����һ����ȵ�k������ӳ٣���ʼ�Ÿɱ�Ϊ5dB�����귽�������ź��Ѷ���
k=3; %�ӳٵ���
r=wgn(1,N,-20); %����-20dB��˹����
xn=wavread('risesun.wav'); %��������������Ŀ���ź��Ѷ��룬�����ź����ӳ�
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
% 
% snr1=0;
% k1=5;
% r=wgn(1,N,-20); %����-20dB��˹����
% xnn=wavread('d:\noisex-92\f16.wav'); %��������������Ŀ���ź��Ѷ��룬�����ź����ӳ�
% xnn=xnn-mean(xnn);
% xnn=xnn/max(abs(xnn));
% xnn=xnn';
% [x1,noise1]=add_noisedata(x1,xnn,fs,fs,snr1); %��һ·��Ӹ��ŵĴ�������
% xnn1=[r(1:k1) xnn(1:N-k1)];
% xnn1=xnn1';
% [x2,noise2]=add_noisedata(x2,xnn1,fs,fs,snr1); %�ڶ�·����ӳٸ��ŵĴ�������
% xnn2=[r(1:2*k1) xnn(1:N-2*k1)];
% xnn2=xnn2';
% [x3,noise3]=add_noisedata(x3,xnn2,fs,fs,snr1); %�ڶ�·����ӳٸ��ŵĴ�������
% xnn3=[r(1:3*k1) xnn(1:N-3*k1)];
% xnn3=xnn3';
% [x4,noise4]=add_noisedata(x4,xnn3,fs,fs,snr1); %�ڶ�·����ӳٸ��ŵĴ�������
% sound(x0,fs);
% sound(xn,fs);
% sound(x1-x2,fs);
% sound(x2-x3,fs);
% sound(x3-x4,fs);
% sound(x4,fs);

% % ������·����ͼ
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
% wlen=256;
% SP=0.5;
% shiftlen=0.5*wlen;
% wnd=hamming(wlen);
% y1=segment(x1,wlen,SP,wnd);    %��֡�����õ�segment����
% y2=segment(x3,wlen,SP,wnd);
% framenum=size(y1,2);           %֡��
% Y1=fft(y1);
% Y2=fft(y2);
% for i=1:framenum
%     for j=1:wlen
%         Pxy(j,i)=Y1(j,i).*conj(Y2(j,i));
%         Pxx(j,i)=Y1(j,i).*conj(Y1(j,i));
%         Pyy(j,i)=Y2(j,i).*conj(Y2(j,i));
%        Txy(j,i)=real(Pxy(j,i)./sqrt(Pxx(j,i).*Pyy(j,i)));
%     end
% end
% 
% T1=0.35;T2=0.8;
% for i=1:framenum
%     for j=1:wlen
%           if Txy(j,i)<T1
%              Txy(j,i)=0.001;
%          elseif Txy(j,i)>T2
%              Txy(j,i)=0.99;
%           end
%     end
% end
% 
% %�õ�����˲���Txy֮�󣬽��зֱ��������˲������ص���ӵõ��µ�x1��x2;
% j=sqrt(-1);
% Y11=Txy.*abs(Y1);
% Y22=Txy.*abs(Y2);
% phase=angle(Y1);
% XXX=exp(j*phase);
% Spec1=Y11.*XXX;
% Spec2=Y22.*XXX;
% y11=zeros((framenum-1)*shiftlen+wlen,1);
% y22=zeros((framenum-1)*shiftlen+wlen,1);
% for i=1:framenum
%     start=(i-1)*shiftlen+1;
%     spec1=Spec1(:,i);
%     spec2=Spec2(:,i);
%     y11(start:start+wlen-1)=y11(start:start+wlen-1)...
%         +real(ifft(spec1,wlen));
%     y22(start:start+wlen-1)=y22(start:start+wlen-1)...
%         +real(ifft(spec2,wlen));
% end
% y3=segment(x2,wlen,SP,wnd);    %��֡�����õ�segment����
% y4=segment(x4,wlen,SP,wnd);
% framenum=size(y1,2);           %֡��
% Y3=fft(y3);
% Y4=fft(y4);
% for i=1:framenum
%     for j=1:wlen
%         Pxy(j,i)=Y1(j,i).*conj(Y2(j,i));
%         Pxx(j,i)=Y1(j,i).*conj(Y1(j,i));
%         Pyy(j,i)=Y2(j,i).*conj(Y2(j,i));
%        Txy(j,i)=real(Pxy(j,i)./sqrt(Pxx(j,i).*Pyy(j,i)));
%     end
% end
% 
% T1=0.35;T2=0.8;
% for i=1:framenum
%     for j=1:wlen
%           if Txy(j,i)<T1
%              Txy(j,i)=0.001;
%          elseif Txy(j,i)>T2
%              Txy(j,i)=0.99;
%           end
%     end
% end
% 
% 
% j=sqrt(-1);
% Y33=Txy.*abs(Y3);
% Y44=Txy.*abs(Y4);
% phase=angle(Y3);
% XXX=exp(j*phase);
% Spec1=Y33.*XXX;
% Spec2=Y44.*XXX;
% y33=zeros((framenum-1)*shiftlen+wlen,1);
% y44=zeros((framenum-1)*shiftlen+wlen,1);
% for i=1:framenum
%     start=(i-1)*shiftlen+1;
%     spec1=Spec1(:,i);
%     spec2=Spec2(:,i);
%     y33(start:start+wlen-1)=y33(start:start+wlen-1)...
%         +real(ifft(spec1,wlen));
%     y44(start:start+wlen-1)=y44(start:start+wlen-1)...
%         +real(ifft(spec2,wlen));
% end
% x1=y11;x2=y22;x3=y33;x4=y44;
% %�����԰����������ǿ�����㷨�������̶������γɺ���������ģ�顣
%d=(y11+y22+y33+y44)/4;
d=(x1+x2+x3+x4)/4;     %�̶������γ�



% frameL=256; % ����֡��
% len=length(s);
% se=zeros(1,len);
% sv=zeros(1,len);
% se_g=zeros(1,len);
% sv_g=zeros(1,len);
% M=33;
% Rx=zeros(M);
% Rs=zeros(M,1);
% x1=x1';x2=x2';x3=x3';x4=x4';d=d';
% for n=1:fix(len/frameL*2-1)
%         Y1=fft([x1((n-1)*frameL/2+1:n*frameL/2+frameL/2).*hanning(frameL)',zeros(1,frameL)]);
%         Y2=fft([x2((n-1)*frameL/2+1:n*frameL/2+frameL/2).*hanning(frameL)',zeros(1,frameL)]);
%         Y3=fft([x3((n-1)*frameL/2+1:n*frameL/2+frameL/2).*hanning(frameL)',zeros(1,frameL)]);
%         Y4=fft([x4((n-1)*frameL/2+1:n*frameL/2+frameL/2).*hanning(frameL)',zeros(1,frameL)]);
%         Ys=fft([d((n-1)*frameL/2+1:n*frameL/2+frameL/2).*hanning(frameL)',zeros(1,frameL)]);
% 
%         A=1/4.*(Y1.*conj(Y1)+Y2.*conj(Y2)+Y3.*conj(Y3)+Y4.*conj(Y4));
%         C=1/6.*(Y1.*conj(Y2)+Y1.*conj(Y3)+Y1.*conj(Y4)+Y2.*conj(Y3)+Y2.*conj(Y4)+Y3.*conj(Y4));
% 
%         W1=C./A;  % �˴���wiener�˲�ԭ�ͼ���
%         SE1=W1.*Ys;
%         se1=real(ifft(SE1));
%         se11=se1(1:frameL);
%         sv((n-1)*frameL/2+1:n*frameL/2+frameL/2)=se11;
%         se=se+sv;
%         sv=zeros(1,len);
%         % ���������µķ���
%         a=real(ifft(A));
%         rr=a(1:M);          
% 
%         for i=1:M
%             for j=1:M
%                 t=abs(i-j);
%                 Rx(i,j)=rr(t+1);
%             end
%         end
%         %Rx=Rx/M;    % Rx Ϊ�����źŵ������ 
%         Cij=real(C);
%         CC=Cij;
%         V=Cij;
%         I=find(Cij<0);
%         Cij(I)=0;
%         J=find(Cij>0);
%         V(J)=0;
%         alpha=Cij.*Cij./(Cij.*Cij+1/6*V);
%         P=alpha.*CC;
%         c=real(ifft(P));
%         Rs=c(1:M)';
%         w=inv(Rx)*Rs;
%         B=0;
%         se1=filter(B,w',d((n-1)*frameL/2+1:n*frameL/2+frameL/2).*hanning(frameL)');
%         sv_g((n-1)*frameL/2+1:n*frameL/2+frameL/2)=se11;
%         se_g=se_g+sv_g;
%         sv_g=zeros(1,len);
%         
% end
% 
% d=se_g;

% �ڹ̶������γ�һ·���ά���˲���
IS=.15; %�޻���ʱ��
output=WienerScalart96m_2(d,fs,IS,0.12);
d=output;
N=length(d);



B=[1,-1,0,0;0,1,-1,0;0,0,1,-1];
% X=[x1';x2';x3';x4'];
% n=B*X;
%����Ӧ��������ģ�飬ʹ�ù�һ����С�����㷨���������������ƣ�

u=0.001;
M=32;
% w=zeros(M,1);
var=zeros(1,M);
%yout=zeros(1,N-1);
w1=zeros(M,1);
w2=zeros(M,1);
w3=zeros(M,1);
y=zeros(1,M);
w=[w1';w2';w3'];
for i=M:N-1
    input1=d(i);
    X1=x1(i:-1:i-M+1);
    X2=x2(i:-1:i-M+1);
    X3=x3(i:-1:i-M+1);
    X4=x4(i:-1:i-M+1);
    X=[X1';X2';X3';X4'];
    Z=B*X;
    Y=w.*Z;
    yout(i)=sum(Y(:));
    e(i)=input1-yout(i);
%     n(i)=sum(Z(:));
%     var(i)=0.2*n(i)^2+0.8*var(i-1);
%     U1=u/var(i);
    w=w+u.*e(i).*Z;
end
% %����Ӧ����ģ�飬�˲�������Ϊ512
% B=[1,-1,0,0;0,1,-1,0;0,0,1,-1]
% u=0.001;
% sysorder=32;
% N=4; 
% M=length(d);
% w1=zeros(sysorder,1);
% w2=zeros(sysorder,1);
% w3=zeros(sysorder,1);
% y=zeros(1,sysorder);
% W=[w1';w2';w3'];
% d=d';
% for n=sysorder+1:M
%     X1=x1(n-sysorder+1:1:n);
%     X2=x2(n-sysorder+1:1:n);
%     X3=x3(n-sysorder+1:1:n);
%     X4=x4(n-sysorder+1:1:n);
%     X=[X1';X2';X3';X4'];
%     Z=B*X;
%     Y=W.*Z;
%     y(n)=sum(Y(:));
%     e(n)=d(n)-y(n);
%     W=W+u*e(n)*Z;
% end
figure(1)
subplot 211;
plot(s);title('Ŀ�������ź�');
sound(s);pause(2);
subplot 212;
plot(xn(1:length(d)));title('�����ź�');
sound(xn(1:length(d)));pause(2);
figure(2)
subplot 211;
plot(x1);
sound(x1);title('��һ·�ɼ����������ź�');pause(2);axis([0 2e4 min(x1) max(x1)]);
subplot 212;
plot(d);title('�̶������γɺ��ź�');
sound(d);pause(2);
figure(3)
subplot 211;
plot(yout);title('���Ƴ����ĸ����ź�');
sound(yout);pause(2);
subplot 212;
plot(e);title('����Ӧ�����γɺ��ź�');
sound(e);pause(2);
wavwrite(e,fs,'none-5');
% m=4;
% Wc=[0.25;0.25;0.25;0.25];
% wop=Wc;
% drawpp(m,wop);
% wop=B'*W(:,sysorder);
% drawpp(m,wop);
% wop=Wc-B'*W(:,sysorder);
% drawpp(m,wop);
% name1=strcat('y1','.pcm');
% fid=fopen(name1,'w');
% fwrite(fid,y(1,:)*32767/3,'int16');
% name2=strcat('out_p','.pcm');
% fid=fopen(name2,'w');
% fwrite(fid,e(1,:)*32767/3,'int16');
% name3=strcat('xn','.pcm');
% fid=fopen(name3,'w');
% fwrite(fid,xn(1,:)*32767/3,'int16');
% x1=x1';
% name4=strcat('x1','.pcm');
% fid=fopen(name4,'w');
% fwrite(fid,x1(1,:)*32767/3,'int16');
% fclose all;
% ����Ϊ��wiener�˲�
e=e';
snr0=SNR_singlech(s(1:length(e)),e)           % ����ά���˲���������
wavwrite(e,fs,'f16-0-0-general')
snr1=SNR_singlech(s(1:length(d)),d)
snr2=SNR_singlech(s(1:length(x1)),x1)
    

% 
% n=y11-y22;         %�����ź�
% %n1=x1-x2;
% %����Ӧ��������ģ�飬ʹ�ù�һ����С�����㷨���������������ƣ�
% 
% u=0.001;
% M=32;
% w=zeros(M,1);
% var=zeros(1,M);
% %yout=zeros(1,N-1);
% for i=M:N-1
%     input1=d(i);
%     input2=n(i:-1:i-M+1);
%     yout(i)=w(1:M)'*input2;
%     e(i)=input1-yout(i);
%     var(i)=0.2*n(i)^2+0.8*var(i-1);
%     U1=u/var(i);
%     w=w+U1.*e(i).*input2;
% end
% %u1=0.0001;var1=zeros(1,M);
% % for i=M:N
% %     input1=d1(i);
% %     input2=n1(i:-1:i-M+1);
% %     yout(i)=w(1:M)*input2;
% %     e1(i)=input1-yout(i);
% %     var1=0.2*(input2)'.^2+0.8*var1;
% %     U1=u./var1;
% %     w=w+U1.*(e1(i).*input2');
% % end
% % 
% 
%  snr1=SNR_singlech(s(1:N-1),d);fprintf(' snr1=%5.1f\n',snr1); 
%  snr2=SNR_singlech(s(1:N-1),e);fprintf(' snr2=%5.1f\n',snr2);
% %  snr3=SNR_singlech(s,d1);fprintf(' snr3=%5.1f\n',snr3); 
% %  snr4=SNR_singlech(s,e1);fprintf(' snr4=%5.1f\n',snr4);
% 
%  sound(y11);
%  pause(1)
%  sound(e);
% %  pause(1);
% %  sound(e1);
% % %  
% figure 
% subplot(511),plot(s);title('clean signal');  axis([1 N -1 1]);
% subplot(512),plot(x1);title('first signal');   axis([1 N -1 1]);
% subplot(513),plot(x2);title('second signal');  axis([1 N -1 1]);
% subplot(514),plot(e);title('enhanced signal'); axis([1 N -1 1]);
% subplot(515),plot(d);title('enhanced signal1'); axis([1 N -1 1]);
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

    

    
    
    
    
