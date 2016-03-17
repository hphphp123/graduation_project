
clc;clear;
load data16k.mat
%设置输入设备参数
% ai = analoginput('mcc');
% addchannel(ai,0:2);
fs=16000;   %采样率
L = 1; %升采样倍数
% set(ai,'samplerate',fs);
T=0.5;      %s: 搜集数据T秒，计算一次位置
% buffer = fs * T; %单次采样时间，对应采样点数
% set(ai,'SamplesPerTrigger',buffer);
% start(ai);
% data = getdata(ai);
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
v=340000;      %mm/s: 音速
NoiseLevel = 0.05;
thta0 = 0; %前一时刻角度值
%麦克风坐标
d=60;    %麦克风间距,mm
Loc_A = [-d,0];
Loc_B = [0,0];
Loc_C = [d,0];
%*用GCC（Generalized Cross-Correlation） Phase Transform 算法估计时延差*
win=kaiser(513,3.4);
fband=[0 200/fs-0.001 200/fs 8000/fs 8000/fs+0.001 1];
m=[0 0 1 1 0 0];
firK = fir2(512,fband,m,win);
parameter_Delay = struct('sample',fs, 'distance',d,'nCal',nMic,'soundspeed',v,'thtapre',thta0,'FIRK',firK);

for ii=1:1
%     start(ai);
    data = data;
    thta000=parameter_Delay.thtapre;
    n=1;
    LL = length(data);
    for i=1:LL
        if (abs(data(i,1))>NoiseLevel)    %%%判断语音涵盖范围，如果语音太少，则不计算%%
            SeNum(n,1) = i;
            n=n+1;
        end
    end
    if(n<LL/10)
        thta = parameter_Delay.thtapre;
    else
       [thta,A,B]= EVD_Delay(data,parameter_Delay);
    end
    parameter_Delay.thtapre=thta;
end

