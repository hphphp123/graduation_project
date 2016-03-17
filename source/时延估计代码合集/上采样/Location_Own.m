
% clc;clear;
%���������豸����
ai = analoginput('mcc');
addchannel(ai,0:2);
fs=16000;   %������
L = 10; %����������
set(ai,'samplerate',fs);
T=0.5;      %s: �Ѽ�����T�룬����һ��λ��
buffer = fs * T; %���β���ʱ�䣬��Ӧ��������
set(ai,'SamplesPerTrigger',buffer);
% start(ai);
% data = getdata(ai);
tMic=0:1/fs:T-1/fs;   %��������ʱ��
nMic=length(tMic);    %�������ݳ���
v=340000;      %mm/s: ����
NoiseLevel = 0.002;
thta0 = 0; %ǰһʱ�̽Ƕ�ֵ
%��˷�����
d=60;    %��˷���,mm
Loc_A = [-d,0];
Loc_B = [0,0];
Loc_C = [d,0];
%*��GCC��Generalized Cross-Correlation�� Phase Transform �㷨����ʱ�Ӳ�*
win=kaiser(513,3.4);
fband=[0 200/fs-0.001 200/fs 8000/fs 8000/fs+0.001 1];
m=[0 0 1 1 0 0];
firK = fir2(512,fband,m,win);

win2=kaiser(513,3.4);
fband2=[0 1/L 1/L+0.001 1];
m2=[1 1 0 0];
firK2 = fir2(512,fband,m,win);
parameter_Delay = struct('sample',fs, 'distance',d,'nCal',nMic,'soundspeed',v,'thtapre',thta0,'FIRK',firK,'Nresample',L,'firKlow',firK2);

for ii=1:1
%     start(ai);
%     data = getdata(ai);
%     data = importdata('badsignal.mat');
    thta000=parameter_Delay.thtapre;
    n=1;
    LL = length(data);
    for i=1:LL
        if (abs(data(i,1))>NoiseLevel)    %%%�ж��������Ƿ�Χ���������̫�٣��򲻼���%%
            SeNum(n,1) = i;
            n=n+1;
        end
    end
    if(n<LL/10)
        thta = parameter_Delay.thtapre;
    else
       [thta,A,B]= GCC_Delay(data,parameter_Delay);
    end
    parameter_Delay.thtapre=thta;
%     thtaout(ii)=thta;
%     Aout(ii)=A;
%     Bout(ii)=B;
end

