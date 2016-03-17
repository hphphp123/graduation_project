
function [thta,A,B] = GCC_Delay(data ,parameter)

nMic = parameter.nCal;
d = parameter.distance;
v = parameter.soundspeed;
fs = parameter.sample;
% thta0 = parameter.thtapre;
thta0 = 360;
thta = parameter.thtapre;
firK = parameter.FIRK;
L = parameter.Nresample;
firKlow = parameter.firKlow;
% Phase Transform�㷨����ʱ��
%����ԭʼ���ݽ�����ʱ����
% L = 1;
% nMic = 1024;
sigMicA_in = data(1:nMic,1);
sigMicB_in = data(1:nMic,2);
sigMicC_in = data(1:nMic,3);
sigMicA = sigMicA_in;
sigMicB = sigMicB_in;
sigMicC = sigMicC_in;

%��ͨ�˲���
sigMicAconv = conv(data(1:nMic,1),firK);
sigMicBconv = conv(data(1:nMic,2),firK);
sigMicCconv = conv(data(1:nMic,3),firK);
% 
sigMicA = sigMicAconv;
sigMicB = sigMicBconv;
sigMicC = sigMicCconv;

%�ز���
% sigMicA_Res = resample(sigMicA_in,L,1,513,3.4);
% sigMicB_Res = resample(sigMicB_in,L,1,513,3.4);
% sigMicC_Res = resample(sigMicC_in,L,1,513,3.4);

if (L>1)
    m = length(sigMicAconv);
    for n=1:(L*(m-1)+1)
    X=(n-1)/L+1;
    if ( X == fix(X) )
        vA(n) = sigMicAconv(X);
        vB(n) = sigMicBconv(X);
        vC(n) = sigMicCconv(X);
    else
        vA(n) = 0;
        vB(n) = 0;
        vC(n) = 0;
    end
    end
    
%     m = length(sigMicA_in);
%     for n=1:(L*(m-1)+1)
%         X=(n-1)/L+1;
%         if ( X == fix(X) )
%            vA(n) = sigMicA_in(X);
%            vB(n) = sigMicB_in(X);
%            vC(n) = sigMicC_in(X);
%         else
%            vA(n) = 0;
%            vB(n) = 0;
%            vC(n) = 0;
%         end
%     end
    sigMicA_Res = conv(vA,firKlow);
    sigMicB_Res = conv(vB,firKlow);
    sigMicC_Res = conv(vC,firKlow);
    
    sigMicA = sigMicA_Res';
    sigMicB = sigMicB_Res';
    sigMicC = sigMicC_Res';
else
    L=1;
end

% GCC_Delay
% FA=fft(sigMicA,length(sigMicA));  %����Ҷ�任
% FB=fft(sigMicB,length(sigMicB));
% FC=fft(sigMicC,length(sigMicC));
% 
% FBC=FB.*conj(FC);               %�󻥹�����
% Amplitude=abs(FB).*abs(FC);     %�����
% GCC=FBC./Amplitude;             %ȥ��������Ϣ
% R_BC=ifft(GCC);
% [val1,DelayDifferBC]=max(R_BC);  %��������ֵ��λ���������ӳٲ�
% 
% FAC=FA.*conj(FC);               %�󻥹�����
% Amplitude=abs(FA).*abs(FC);     %�����
% GCC=FAC./Amplitude;             %ȥ��������Ϣ
% R_AC=ifft(GCC);
% [val1,DelayDifferAC]=max(R_AC);  %��������ֵ��λ���������ӳٲ
% DelayDifferAC = DelayDifferAC
% 
% FAB=FA.*conj(FB);               %�󻥹�����
% Amplitude=abs(FA).*abs(FB);     %�����
% GCC=FAB./Amplitude;             %ȥ��������Ϣ
% R_AB=ifft(GCC);
% [val1,DelayDifferAB]=max(R_AB);  %��������ֵ��λ���������ӳٲ
% DelayDifferAB = DelayDifferAB
% 
% datanBC = fs*L*d/v*2;
% datanAC = fs*L*d/v*2*2;
% if(DelayDifferBC<datanBC)
%     distDiffBC=(DelayDifferBC-1)/fs/L*v;
% else distDiffBC=(DelayDifferBC-length(sigMicA)-1)/fs/L*v;
% end
% 
% if(DelayDifferAC<datanAC)
%     distDiffAC=(DelayDifferAC-1)/fs/L*v;
% else distDiffAC=(DelayDifferAC-length(sigMicA)-1)/fs/L*v;
% end
% 
% if(DelayDifferAB<datanBC)
%     distDiffAB=(DelayDifferAB-1)/fs/L*v;
% else distDiffAB=(DelayDifferAB-length(sigMicA)-1)/fs/L*v;
% end
% 
% A=distDiffBC;
% B=distDiffAC;
% C=distDiffAB;

%%NLMS
x1=sigMicA ;
x2=sigMicB;
x3=sigMicC;

len=length(x1);
x=[zeros(1000,1);x1];
w=zeros(1000,1);%��ʼ1000�׼�Ȩϵ�����˲���������Ϊ1000
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
else delay=(delay-length(sigMicA));
end
delay_lms1=delay/fs/L*v;

len=length(x2);
x=[zeros(1000,1);x2];
w=zeros(1000,1);%��ʼ1000�׼�Ȩϵ�����˲���������Ϊ1000
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
else delay=(delay-length(sigMicA));
end
delay_lms2=delay/fs/L*v;

A=delay_lms2
B=delay_lms1

%%�ж��Ƿ���0�Ⱥ�180�ȣ���ʱ���㹫ʽ������%%%
if(B ==(2*A)&&A>0)
    thtaround=abs(thta0-0);
    dround=abs(abs(B)-2*d);
    if(thtaround<10 && dround<5)
        thta = 0;
    else
        B=B+0.01;
        %%%����A��B��˴�������r3/r2%%%
        M=2*A^2+B^2-4*A*B+2*d^2;
        N=-2*(2*A-B);
        r3=M/N;
        r2=A+r3-B;
        cosA=(r3^2+d^2-r2^2)/(2*d*r3);
        cosA = min(cosA,1);
        thtaCal = acos(cosA)*180/pi;
        if(r3>0&&r2>0)
            thta = thtaCal;
        end
    end
elseif(B == (2*A)&&A<0)
    thtaround=abs(thta0-180);
    dround=abs(abs(B)-2*d);
    if(thtaround<10&& dround<5)
        thta = 180;
    else
        B=B+0.01;
        %%%����A��B��˴�������r3/r2%%%
        M=2*A^2+B^2-4*A*B+2*d^2;
        N=-2*(2*A-B);
        r3=M/N;
        r2=A+r3-B;
        cosA=(r3^2+d^2-r2^2)/(2*d*r3);
        cosA = min(cosA,1);
        thtaCal = acos(cosA)*180/pi;
        if(r3>0&&r2>0)
            thta = thtaCal;
        end 
    end
else
    %%%����A��B��˴�������r3/r2%%%
    M=2*A^2+B^2-4*A*B+2*d^2;
    N=-2*(2*A-B);
    r3=M/N;
    r2=A+r3-B;
    cosA=(r3^2+d^2-r2^2)/(2*d*r3);
    cosA = min(cosA,1);
    thtaCal = acos(cosA)*180/pi;
    if(r3>0&&r2>0)
        thta = thtaCal;
    else
        Around = abs(abs(A)-d);
        Bround = abs(abs(B)-2*d);
        if(Around<10&&Bround<20&&A<0&&B<0)
            thta =180;
        elseif(Around<20&&Bround<10&&A<0&&B<0)
            thta =180;
        elseif(Around<10&&Bround<20&&A>0&&B>0)
            thta =0;
        else(Around<20&&Bround<10&&A>0&&B>0)
            thta =0;
        end
    end
end
%     sinA=sqrt(1-cosA^2);
    %
%     a=Loc_A(1)+r3*cosA;
%     b=r3*sinA;
%     Loc_reslut=[a,b];
