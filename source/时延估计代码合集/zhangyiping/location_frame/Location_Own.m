
clc;clear;
%���������豸����
ai = analoginput('mcc');
addchannel(ai,0:2);
fs=100000;   %������
L = 1; %����������
set(ai,'samplerate',fs);
T=2.5;      %s: �Ѽ�����T�룬����һ��λ��
buffer = fs * T; %���β���ʱ�䣬��Ӧ��������
set(ai,'SamplesPerTrigger',buffer);

tMic=0:1/fs:T-1/fs;   %��������ʱ��
nMic=length(tMic);    %�������ݳ���
Lframe = 512;
Nframe = fix(nMic/Lframe);

v=340000;      %mm/s: ����
NoiseLevel = 0.0002;
thta0 = 360; %ǰһʱ�̽Ƕ�ֵ
%��˷�����
d=60;    %��˷���,mm
Loc_A = [-d,0];
Loc_B = [0,0];
Loc_C = [d,0];
%*��GCC��Generalized Cross-Correlation�� Phase Transform �㷨����ʱ�Ӳ�*
win=kaiser(513,3.4);
fband=[0 50/fs-0.0001 50/fs 8000/fs 8000/fs+0.0001 1];
m=[0 0 1 1 0 0];
firK = fir2(512,fband,m,win);

win2=kaiser(513,3.4);
fband2=[0 1/L 1/L+0.001 1];
m2=[1 1 0 0];
firK2 = fir2(512,fband,m,win);
parameter_Delay = struct('sample',fs, 'distance',d,'nCal',Lframe,'soundspeed',v,'thtapre', thta0,...
                        'FIRK',firK,'Nresample',L,'firKlow',firK2);

stop=0;
ii = 1;
while(stop==0)
    if(ii>1)
        stop=1;
    else
        start(ai);
        datain = getdata(ai);
%         dclock = fix(clock);
%         nameclock=strcat(num2str(dclock(1,1)),'year',num2str(dclock(1,2)),num2str(dclock(1,3)),'month',num2str(dclock(1,4)),num2str(dclock(1,5)));
%         name1=strcat('channel1_',nameclock,'.pcm');
%         name2=strcat('channel2_',nameclock,'.pcm');
%         name3=strcat('channel3_',nameclock,'.pcm');
%         fid = fopen(name1,'w');
%         fwrite(fid,datain(:,1)*32767,'int16');
%         fid = fopen(name2,'w');
%         fwrite(fid,datain(:,2)*32767,'int16');
%         fid = fopen(name3,'w');
%         fwrite(fid,datain(:,3)*32767,'int16');
        % datain = importdata('badsignal.mat');
        thta000=parameter_Delay.thtapre;
        for jj=1:Nframe
            ncut = (jj-1)*Lframe+1;
            data(:,1) = datain(ncut:(ncut+Lframe-1),1);
            data(:,2) = datain(ncut:(ncut+Lframe-1),2);
            data(:,3) = datain(ncut:(ncut+Lframe-1),3);
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
                [thta,A,B]= LMS_Delay(data,parameter_Delay);
            end
            
            thtaout(jj,ii)=thta;
            Aout(jj,ii)=A;
            Bout(jj,ii)=B;    
        end
        parameter_Delay.thtapre=thta;
        ii=ii+1;
    end
end

