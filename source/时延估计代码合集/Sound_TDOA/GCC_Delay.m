
function [thta,A,B] = GCC_Delay(data ,parameter)

nMic = parameter.nCal;
d = parameter.distance;
v = parameter.soundspeed;
fs = parameter.sample;
thta0 = parameter.thtapre;
thta = parameter.thtapre;
firK = parameter.FIRK;
% Phase Transform�㷨����ʱ��
sigMicA_in = data(1:nMic,1);
sigMicB_in = data(1:nMic,2);
sigMicC_in = data(1:nMic,3);
%��ͨ�˲�
sigMicAconv = conv(data(1:nMic,1),firK);
sigMicBconv = conv(data(1:nMic,2),firK);
sigMicCconv = conv(data(1:nMic,3),firK);

sigMicA = sigMicAconv;
sigMicB = sigMicBconv;
sigMicC = sigMicCconv;

% sigMicA_Con = sigMicAconv(513:nMic+512,1);
% sigMicB_Con = sigMicBconv(513:nMic+512,1);
% sigMicC_Con = sigMicCconv(513:nMic+512,1);

%�ز���
L=1;
% sigMicA_Res = resample(sigMicA_in,L,1,513,3.4);
% sigMicB_Res = resample(sigMicB_in,L,1,513,3.4);
% sigMicC_Res = resample(sigMicC_in,L,1,513,3.4);

% sigMicA_inter = interp(sigMicA_in,L);
% sigMicB_inter = interp(sigMicB_in,L);
% sigMicC_inter = interp(sigMicC_in,L);

%     sigMicA = sigMicA_in;
%     sigMicB = sigMicB_in;
%     sigMicC = sigMicC_in;

% sigMicA = sigMicA_in;
% sigMicB = sigMicB_in;
% sigMicC = sigMicC_in;
% sig(:,1) = sigMicA;
% sig(:,2) = sigMicB;
% sig(:,3) = sigMicC;

% FAin=abs(fft(sigMicAin,length(sigMicAin)));  %����Ҷ�任
% FAA=abs(fft(sigMicA,length(sigMicA))); 
% FCin=fft(sigMicCin,length(sigMicCin));

FA=fft(sigMicA,length(sigMicA));  %����Ҷ�任
FB=fft(sigMicB,length(sigMicB));
FC=fft(sigMicC,length(sigMicC));

FBC=FB.*conj(FC);               %�󻥹�����
Amplitude=abs(FB).*abs(FC);     %�����
GCC=FBC./Amplitude;             %ȥ��������Ϣ
R_BC=ifft(GCC);
[val1,DelayDifferBC]=max(R_BC);  %��������ֵ��λ���������ӳٲ�

FAC=FA.*conj(FC);               %�󻥹�����
Amplitude=abs(FA).*abs(FC);     %�����
GCC=FAC./Amplitude;             %ȥ��������Ϣ
R_AC=ifft(GCC);
[val1,DelayDifferAC]=max(R_AC);  %��������ֵ��λ���������ӳٲ

FAB=FA.*conj(FB);               %�󻥹�����
Amplitude=abs(FA).*abs(FB);     %�����
GCC=FAB./Amplitude;             %ȥ��������Ϣ
R_AB=ifft(GCC);
[val1,DelayDifferAB]=max(R_AB);  %��������ֵ��λ���������ӳٲ

datanBC = fs*L*d/v*3;
datanAC = fs*L*d/v*3*2;
if(DelayDifferBC<datanBC)
    distDiffBC=(DelayDifferBC-1)/fs*v;
else distDiffBC=(DelayDifferBC-length(sigMicA)*L-1)/fs*v;
end

if(DelayDifferAC<datanAC)
    distDiffAC=(DelayDifferAC-1)/fs*v;
else distDiffAC=(DelayDifferAC-length(sigMicA)*L-1)/fs*v;
end

if(DelayDifferAB<datanBC)
    distDiffAB=(DelayDifferAB-1)/fs*v;
else distDiffAB=(DelayDifferAB-length(sigMicA)*L-1)/fs*v;
end

A=distDiffBC;
B=distDiffAC;
C=distDiffAB;

%%�ж��Ƿ���0�Ⱥ�180�ȣ���ʱ���㹫ʽ������%%%
if(B == (2*A)&&A>0)
    thtaround=abs(thta0-0);
    if(thtaround<10)
        thta = 0;
    end
elseif(B == (2*A)&&A<0)
    thtaround=abs(thta0-180);
    if(thtaround<10)
        thta = 180;
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
    end
end
%     sinA=sqrt(1-cosA^2);
    %
%     a=Loc_A(1)+r3*cosA;
%     b=r3*sinA;
%     Loc_reslut=[a,b];
