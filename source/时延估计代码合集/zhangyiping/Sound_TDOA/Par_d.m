clc
clear all
close all
err_data=[];
Error_Avr_data=[];
%%
% *����������*
for d = 10:40
x=5:5:100;
y=20:5:80;
[xx,yy]=meshgrid(x, y); % xx��yy����25x25�ľ���
zz=xx.*0; % ���㺯��ֵ��zzҲ��21x21�ľ���  
for  x1=5:5:100 
    for y1=20:5:80
%--��Դ��ز���
fmin=500;
fmax=2000;    %Hz: ��ԴΪһƵ�ʽ���������źţ����Ƶ��fmin�����Ƶ��fmax
S_last=0.1;   %s ����Դ����ʱ��
%--�������źŴ�����ز���
fs=3e6;   %������ 
ts=1/fs;    %�������

T=0.12;      %s: �Ѽ�����T�룬����һ��λ��
tMic=0:1/fs:T-1/fs;   %��������ʱ��
nMic=length(tMic);    %�������ݳ���
%--�����������
v=34000;      %cm/s: ����
%--��Դ��MICλ��   ��λCM
% d=20;    %��˷���
Lco_S = [x1,y1];
Loc_A = [50-d,0];
Loc_B = [50,0];
Loc_C = [50+d,0];
%%
% *ģ����Դ�źż���MIC���յ����ź�*

%--������Դ
t = 0:ts:S_last; %������������ʱ��0.2�� 
s = chirp(t,fmin,S_last,fmax,'linear');%Դ�ź�, Ƶ�����Ե����������ź�
nsource=length(s);

%--������Դ������MIC֮����ӳٵ���
diff_SA = sqrt(sum((Lco_S-Loc_A).^2))/v/ts;
diff_SB = sqrt(sum((Lco_S-Loc_B).^2))/v/ts;
diff_SC = sqrt(sum((Lco_S-Loc_C).^2))/v/ts;

%--MIC���յ�������
sigMicA=zeros(1,nMic);
sigMicB=zeros(1,nMic);
sigMicC=zeros(1,nMic);

sigMicA(1+diff_SA:nsource+diff_SA)=s;
sigMicB(1+diff_SB:nsource+diff_SB)=s;
sigMicC(1+diff_SC:nsource+diff_SC)=s;

%%
% *��GCC��Generalized Cross-Correlation�� Phase Transform �㷨����ʱ�Ӳ�*

% Phase Transform�㷨����ʱ��
FB=fft(sigMicB,length(sigMicB));  %����Ҷ�任
FC=fft(sigMicC,length(sigMicC));
FBC=FB.*conj(FC);               %�󻥹�����
Amplitude=abs(FB).*abs(FC);     %�����
GCC=FBC./Amplitude;             %ȥ��������Ϣ
R_BC=ifft(GCC);                 
[val1,DelayDifferBC]=max(R_BC);  %��������ֵ��λ���������ӳٲ

FA=fft(sigMicA,length(sigMicA));  %����Ҷ�任
FC=fft(sigMicC,length(sigMicC));
FAC=FA.*conj(FC);               %�󻥹�����
Amplitude=abs(FA).*abs(FC);     %�����
GCC=FAC./Amplitude;             %ȥ��������Ϣ
R_AC=ifft(GCC);                 
[val1,DelayDifferAC]=max(R_AC);  %��������ֵ��λ���������ӳٲ

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
% %--��Դλ��ͼ
% figure();
% plot(Lco_S(1),Lco_S(2),'ro');
% hold on;
% plot(Loc_reslut(1),Loc_reslut(2),'o');
%--�������
Error=sqrt(sum((Loc_reslut-Lco_S).^2)); %�������

zz(fix((y1-15)/5),fix(x1/5))=Error;
err_data=[err_data Error];
    end
end
mesh(xx, yy, zz); % ����������״ͼ
xlabel('����X'),ylabel('����Y'),zlabel('���ֵ');
Error_Avr = mean(err_data);
Error_Avr_data=[Error_Avr_data Error_Avr];
end
       figure();
       plot(Error_Avr_data,'-*');
       ylabel('���')
       xlabel('d')
       grid on;
       box on;        
set(gcf,'Units','centimeters','Position',[10 10 6.67 5]);%
set(gcf,'Units','centimeters','Position',[10 10 9 6.75]);%
%set(gcf,'Units','centimeters','Position',[10 10 13.5 9]);%

set(get(gca,'XLabel'),'FontSize',8);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',8);
set(get(gca,'TITLE'),'FontSize',8);
set(gca,'fontsize',8);
set(gca,'linewidth',0.5); %�����ߴ�0.5��
%set(gca,'box','off');%Controls the box around the plotting area
set(get(gca,'Children'),'linewidth',1.2);%����ͼ���߿�1.5��















