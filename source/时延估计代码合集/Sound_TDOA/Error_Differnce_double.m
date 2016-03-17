clc
clear all
close all
err_data=[];
%%
% *各参数设置*
for  x=0:1:100 
    
%--声源相关参数
fmin=500;
fmax=2000;    %Hz: 信源为一频率渐变的余弦信号，最低频率fmin，最高频率fmax
S_last=0.2;   %s ：声源持续时间
%--采样和信号处理相关参数
fs=3e6;   %采样率 
ts=1/fs;    %采样间隔

T=0.3;      %s: 搜集数据T秒，计算一次位置
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
%--物理参数设置
v=34000;      %cm/s: 音速
%--声源和MIC位置   单位CM
d=20;    %麦克风间距
Lco_S = [50,x];
Loc_A = [0,0];
Loc_B = [d,0];
Loc_C = [2*d,0];

Loc_E = [100,100];
Loc_F = [100-d,100];
Loc_G = [100-2*d,100];
%%
% *模拟声源信号及各MIC接收到的信号*

%--产生声源
t = 0:ts:S_last; %假设声波持续时间0.2秒 
s = chirp(t,fmin,S_last,fmax,'linear');%源信号, 频率线性递增的余弦信号
nsource=length(s);

%--计算信源传到各MIC之间的延迟点数
diff_SA = sqrt(sum((Lco_S-Loc_A).^2))/v/ts;
diff_SB = sqrt(sum((Lco_S-Loc_B).^2))/v/ts;
diff_SC = sqrt(sum((Lco_S-Loc_C).^2))/v/ts;

diff_SE = sqrt(sum((Lco_S-Loc_E).^2))/v/ts;
diff_SF = sqrt(sum((Lco_S-Loc_F).^2))/v/ts;
diff_SG = sqrt(sum((Lco_S-Loc_G).^2))/v/ts;

%--MIC接收到的数据
sigMicA=zeros(1,nMic);
sigMicB=zeros(1,nMic);
sigMicC=zeros(1,nMic);

sigMicA(1+diff_SA:nsource+diff_SA)=s;
sigMicB(1+diff_SB:nsource+diff_SB)=s;
sigMicC(1+diff_SC:nsource+diff_SC)=s;

sigMicE=zeros(1,nMic);
sigMicF=zeros(1,nMic);
sigMicG=zeros(1,nMic);

sigMicE(1+diff_SE:nsource+diff_SE)=s;
sigMicF(1+diff_SF:nsource+diff_SF)=s;
sigMicG(1+diff_SG:nsource+diff_SG)=s;

%%
% *用GCC（Generalized Cross-Correlation） Phase Transform 算法估计时延差*   

% 第一组传感器

% Phase Transform算法求延时差
FA=fft(sigMicA,length(sigMicA));  %傅里叶变换
FB=fft(sigMicB,length(sigMicB));
FAB=FA.*conj(FB);               %求互功率谱
Amplitude=abs(FA).*abs(FB);     %求幅度
GCC=FAB./Amplitude;             %去除幅度信息
R_ab=ifft(GCC);                 
[val1,DelayDifferAB]=max(R_ab);  %互相关最大值的位置体现了延迟差。

FC=fft(sigMicC,length(sigMicC));
FAC=FA.*conj(FC);               %求互功率谱
Amplitude=abs(FA).*abs(FC);     %求幅度
GCC=FAC./Amplitude;             %去除幅度信息
R_ac=ifft(GCC);                 
[val2,DelayDifferAC]=max(R_ac);  %互相关最大值的位置体现了延迟差。

if(DelayDifferAB<3000)
  distDiffAB=(DelayDifferAB-1)/fs*v;
else distDiffAB=(DelayDifferAB-nMic-1)/fs*v;
end

if(DelayDifferAC<3000)
  distDiffAC=(DelayDifferAC-1)/fs*v;
else distDiffAC=(DelayDifferAC-nMic-1)/fs*v;
end

flag=0;
for a=10:20:100
    for b=10:20:100
%方程求解
funSrc=inline('[sqrt(sum((x-A).^2))-sqrt(sum((x-B).^2))-dAB;sqrt(sum((x-A).^2))-sqrt(sum((x-C).^2))-dAC]','x','A','B','C','dAB','dAC');
[Loc_reslut1,err]=fsolve(funSrc,[a,b],[100,100],Loc_A,Loc_B,Loc_C,distDiffAB,distDiffAC);    %迭代初值为其中一个MIC的位置
Error1=sqrt(sum((Loc_reslut1-Lco_S).^2)); %估计误差
 if Error1<1;
    flag=1;
    break;    
 end
    end
    if flag ==1;
      break;
    end
end

%%
% *用GGG（Generalized Gross-Gorrelation） Phase Transform 算法估计时延差*

% 第二组传感器

% Phase Transform算法求延时差
FE=fft(sigMicE,length(sigMicE));  %傅里叶变换
FF=fft(sigMicF,length(sigMicF));
FEF=FE.*conj(FF);               %求互功率谱
amplitude=abs(FE).*abs(FF);     %求幅度
GCC=FEF./amplitude;             %去除幅度信息
R_ab=ifft(GCC);                 
[val1,DelayDifferEF]=max(R_ab);  %互相关最大值的位置体现了延迟差。

FG=fft(sigMicG,length(sigMicG));
FEG=FE.*conj(FG);               %求互功率谱
amplitude=abs(FE).*abs(FG);     %求幅度
GCC=FEG./amplitude;             %去除幅度信息
R_ac=ifft(GCC);                 
[val2,DelayDifferEG]=max(R_ac);  %互相关最大值的位置体现了延迟差。

if(DelayDifferEF<3000)
  distDiffEF=(DelayDifferEF-1)/fs*v;
else distDiffEF=(DelayDifferEF-nMic-1)/fs*v;
end

if(DelayDifferEG<3000)
  distDiffEG=(DelayDifferEG-1)/fs*v;
else distDiffEG=(DelayDifferEG-nMic-1)/fs*v;
end

flag=0;
for a=10:20:100
    for b=10:20:100
%方程求解
funSrc=inline('[sqrt(sum((x-E).^2))-sqrt(sum((x-F).^2))-dEF;sqrt(sum((x-E).^2))-sqrt(sum((x-G).^2))-dEG]','x','E','F','G','dEF','dEG');
[Loc_reslut2,err]=fsolve(funSrc,[a,b],[100,100],Loc_E,Loc_F,Loc_G,distDiffEF,distDiffEG);    %迭代初值为其中一个MIG的位置
Error2=sqrt(sum((Loc_reslut2-Lco_S).^2)); %估计误差
 if Error2<1;
    flag=1;
    break;    
 end
    end
    if flag ==1;
      break;
    end
end
%%
% *取舍正确的距离*

if ( Error1<1 &&  Error2<1)
   Loc_reslut=(Loc_reslut1+Loc_reslut2)/2;
else if (Error2<1)
    Loc_reslut=Loc_reslut2;
     else Loc_reslut=Loc_reslut1;
    end
end 

%--计算误差
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %估计误差
err_data=[err_data Error];

end
       figure();
       plot(0:100,err_data,'-*');
       ylabel('误差')
       xlabel('横坐标50下的纵坐标')
       grid on;
       box on;
       
          
set(gcf,'Units','centimeters','Position',[10 10 6.67 5]);%
set(gcf,'Units','centimeters','Position',[10 10 9 6.75]);%
%set(gcf,'Units','centimeters','Position',[10 10 13.5 9]);%

set(get(gca,'XLabel'),'FontSize',8);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',8);
set(get(gca,'TITLE'),'FontSize',8);
set(gca,'fontsize',8);
set(gca,'linewidth',0.5); %坐标线粗0.5磅
%set(gca,'box','off');%Controls the box around the plotting area
set(get(gca,'Children'),'linewidth',1.2);%设置图中线宽1.5磅















