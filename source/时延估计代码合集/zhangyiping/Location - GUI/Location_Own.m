
function thta = Location_Own(fs,T,ai,thta0)
%% 定位主函数
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
v=340000;      %mm/s: 音速
NoiseLevel = 0.01;
%麦克风坐标
d=60;    %麦克风间距,mm
% Loc_A = [-d,0];
% Loc_B = [0,0];
% Loc_C = [d,0];
%*用GCC（Generalized Cross-Correlation） Phase Transform 算法估计时延差*
win=kaiser(513,3.4);
fband=[0 50/fs-0.0001 50/fs 8000/fs 8000/fs+0.0001 1];
m=[0 0 1 1 0 0];
firK = fir2(512,fband,m,win);

parameter_Delay = struct('sample',fs, 'distance',d,'nCal',nMic,'soundspeed',v,'thtapre', thta0,...
                        'FIRK',firK);
%% 采集卡数据读取
start(ai);
datain = getdata(ai);
datain(:,1) = datain(:,1)-mean(datain(:,1));
datain(:,2) = datain(:,2)-mean(datain(:,2));
datain(:,3) = datain(:,3)-mean(datain(:,3));
data = datain;
%         data(:,1) = noise_estimation(datain(:,1),fs);
%         data(:,2) = noise_estimation(datain(:,2),fs);
%         data(:,3) = noise_estimation(datain(:,3),fs);
%         dclock = fix(clock);
%         nameclock=strcat(num2str(dclock(1,1)),'year',num2str(dclock(1,2)),num2str(dclock(1,3)),'month',num2str(dclock(1,4)),num2str(dclock(1,5)));
%         name1=strcat('bg8kchannel1_',nameclock,'.pcm');
%         name2=strcat('bg8kchannel2_',nameclock,'.pcm');
%         name3=strcat('bg8kchannel3_',nameclock,'.pcm');
%         fid = fopen(name1,'w');
%         fwrite(fid,datain(:,1)*32767,'int16');    
%         fid = fopen(name2,'w');
%         fwrite(fid,datain(:,2)*32767,'int16');
%         fid = fopen(name3,'w');
%         fwrite(fid,datain(:,3)*32767,'int16');
%         fclose all;
        % datain = importdata('badsignal.mat');
        
n=1;
LL = length(data);
for i=1:LL
    if (abs(data(i,1))>NoiseLevel)    %%%判断语音涵盖范围，如果语音太少，则不计算%%
        n=n+1;
    end
end
if(n<LL/10)
    thta = parameter_Delay.thtapre;
else
    thta= GCC_Delay(data,parameter_Delay);
end

% if(parameter_Delay.thtapre~=360 && abs(thta-parameter_Delay.thtapre)>50)
%     thta = parameter_Delay.thtapre;
% end
parameter_Delay.thtapre=thta;
thtaout = thta



