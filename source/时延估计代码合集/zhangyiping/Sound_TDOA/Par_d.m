clc
clear all
close all
err_data=[];
Error_Avr_data=[];
%%
% *各参数设置*
for d = 10:40
x=5:5:100;
y=20:5:80;
[xx,yy]=meshgrid(x, y); % xx和yy都是25x25的矩阵
zz=xx.*0; % 计算函数值，zz也是21x21的矩阵  
for  x1=5:5:100 
    for y1=20:5:80
%--声源相关参数
fmin=500;
fmax=2000;    %Hz: 信源为一频率渐变的余弦信号，最低频率fmin，最高频率fmax
S_last=0.1;   %s ：声源持续时间
%--采样和信号处理相关参数
fs=3e6;   %采样率 
ts=1/fs;    %采样间隔

T=0.12;      %s: 搜集数据T秒，计算一次位置
tMic=0:1/fs:T-1/fs;   %接收数据时间
nMic=length(tMic);    %接收数据长度
%--物理参数设置
v=34000;      %cm/s: 音速
%--声源和MIC位置   单位CM
% d=20;    %麦克风间距
Lco_S = [x1,y1];
Loc_A = [50-d,0];
Loc_B = [50,0];
Loc_C = [50+d,0];
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

%--MIC接收到的数据
sigMicA=zeros(1,nMic);
sigMicB=zeros(1,nMic);
sigMicC=zeros(1,nMic);

sigMicA(1+diff_SA:nsource+diff_SA)=s;
sigMicB(1+diff_SB:nsource+diff_SB)=s;
sigMicC(1+diff_SC:nsource+diff_SC)=s;

%%
% *用GCC（Generalized Cross-Correlation） Phase Transform 算法估计时延差*

% Phase Transform算法求延时差
FB=fft(sigMicB,length(sigMicB));  %傅里叶变换
FC=fft(sigMicC,length(sigMicC));
FBC=FB.*conj(FC);               %求互功率谱
Amplitude=abs(FB).*abs(FC);     %求幅度
GCC=FBC./Amplitude;             %去除幅度信息
R_BC=ifft(GCC);                 
[val1,DelayDifferBC]=max(R_BC);  %互相关最大值的位置体现了延迟差。

FA=fft(sigMicA,length(sigMicA));  %傅里叶变换
FC=fft(sigMicC,length(sigMicC));
FAC=FA.*conj(FC);               %求互功率谱
Amplitude=abs(FA).*abs(FC);     %求幅度
GCC=FAC./Amplitude;             %去除幅度信息
R_AC=ifft(GCC);                 
[val1,DelayDifferAC]=max(R_AC);  %互相关最大值的位置体现了延迟差。

if(DelayDifferBC<5000)
  distDiffBC=(DelayDifferBC-1)/fs*v;
else distDiffBC=(DelayDifferBC-nMic-1)/fs*v;
end

if(DelayDifferAC<5000)
  distDiffAC=(DelayDifferAC-1)/fs*v;
else distDiffAC=(DelayDifferAC-nMic-1)/fs*v;
end

A=distDiffBC;
B=distDiffAC;

M=4*A^3-B^3+6*A*B^2-10*A^2*B+4*d^2*A-2*d^2*B;
N=-8*A^2-2*B^2+8*A*B;
r3=M/N;

M1=4*A^3-B^3+4*A*B^2-6*A^2*B+2*d^2*B-4*d^2*A;
N1=8*A^2+2*B^2-8*A*B;
r2=M1/N1;

cos=(r3^2+d^2-r2^2)/(2*d*r3);
sin=sqrt(1-cos^2);

a=(50-d)+r3*cos;
b=r3*sin;
Loc_reslut=[a,b];
% %--声源位置图
% figure();
% plot(Lco_S(1),Lco_S(2),'ro');
% hold on;
% plot(Loc_reslut(1),Loc_reslut(2),'o');
%--计算误差
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %估计误差

zz(fix((y1-15)/5),fix(x1/5))=Error;
err_data=[err_data Error];
    end
end
mesh(xx, yy, zz); % 画出立体网状图
xlabel('坐标X'),ylabel('坐标Y'),zlabel('误差值');
Error_Avr = mean(err_data);
Error_Avr_data=[Error_Avr_data Error_Avr];
end
       figure();
       plot(Error_Avr_data,'-*');
       ylabel('误差')
       xlabel('d')
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















