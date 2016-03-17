
function [thta,A,B] = GCC_Delay(data ,parameter)

nMic = parameter.nCal;
d = parameter.distance;
v = parameter.soundspeed;
fs = parameter.sample;
thta0 = parameter.thtapre;
thta = parameter.thtapre;
firK = parameter.FIRK;
% Phase Transform算法求延时差
sigMicA_in = data(1:nMic,1);
sigMicB_in = data(1:nMic,2);
sigMicC_in = data(1:nMic,3);
%带通滤波
sigMicAconv = conv(data(1:nMic,1),firK);
sigMicBconv = conv(data(1:nMic,2),firK);
sigMicCconv = conv(data(1:nMic,3),firK);

sigMicA = sigMicAconv;
sigMicB = sigMicBconv;
sigMicC = sigMicCconv;

% sigMicA_Con = sigMicAconv(513:nMic+512,1);
% sigMicB_Con = sigMicBconv(513:nMic+512,1);
% sigMicC_Con = sigMicCconv(513:nMic+512,1);

%重采样
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

% FAin=abs(fft(sigMicAin,length(sigMicAin)));  %傅里叶变换
% FAA=abs(fft(sigMicA,length(sigMicA))); 
% FCin=fft(sigMicCin,length(sigMicCin));



%%% 分帧+内插，求时延估计delay的众数
x1=sigMicA ;
x2=sigMicB;
x3=sigMicC;

N=2048;
Lov = 4;
M = round(N/Lov);    
%加窗
a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1)));

Nx=length(x1);
Ndata=1+ceil((Nx-N)/M); %帧数
Nf=N;  %傅里叶变换长度
Nfh=Nf/2+1;   %频谱共轭对称，取一半
Sxy=zeros(Nf,1);  %互功率谱
Sx=zeros(Nf,1);  %功率谱
Sy=zeros(Nf,1);
delay=zeros(Ndata,1);  %每帧延迟矩阵
d=abs(d); %距离为d
Nd = 2+ceil(d/v*fs);%最大延迟点数

alpha = 0.8;    %遗忘概率，语音信号非平稳                          
alpha1 = 1-alpha;

m = 0;
for n = 1:M:Nx-N+1
  m = m+1;
  n1 = n:n+N-1; 
  X2 = fft(x2(n1).*w,Nf);
  X3 = fft(x3(n1).*w,Nf);
  X2 = X2(1:Nf);                          
  X3 = X3(1:Nf);
  Sxy =alpha*Sxy + alpha1* X2 .* conj(X3); 
  Sx =alpha*Sx + alpha1* X2 .* conj(X2);
  Sy =  alpha*Sy + alpha1*X3 .* conj(X3);
  Cxy = ifft(Sxy./(abs(Sx).*abs(Sy)),Nf);  
  Cxy = [Cxy(Nf-Nd+1:Nf) ; Cxy(1:Nd)];   %只能在前和后找，要不然就不对了
  [Cxymax,imax] = max(Cxy);
  if imax<Nd
    delay(m) = Nd-imax+1;
  else 
    delay(m)=Nd*2-imax+1;
  end
end
DelayDifferBC=mode(delay);

m = 0;
for n = 1:M:Nx-N+1
  m = m+1;
  n1 = n:n+N-1; 
  X1 = fft(x1(n1).*w,Nf);
  X3 = fft(x3(n1).*w,Nf);
  X1 = X1(1:Nf);                          
  X3 = X3(1:Nf);
  Sxy =alpha*Sxy + alpha1* X1 .* conj(X3); 
  Sx =alpha*Sx + alpha1* X1 .* conj(X1);
  Sy =  alpha*Sy + alpha1*X3 .* conj(X3);
  Cxy = ifft(Sxy./(abs(Sx).*abs(Sy)),Nf);  
  Cxy = [Cxy(Nf-Nd+1:Nf) ; Cxy(1:Nd)];   %只能在前和后找，要不然就不对了
  [Cxymax,imax] = max(Cxy);
  if imax<Nd
    delay(m) = Nd-imax+1;
  else 
    delay(m)=Nd*2-imax+1;
  end
end
DelayDifferAC=mode(delay);


datanBC = fs*L*d/v*3;
datanAC = fs*L*d/v*3*2;
if(DelayDifferBC<datanBC)
    distDiffBC=(DelayDifferBC)/fs*v;
else distDiffBC=(DelayDifferBC-length(sigMicA)*L)/fs*v;
end

if(DelayDifferAC<datanAC)
    distDiffAC=(DelayDifferAC)/fs*v;
else distDiffAC=(DelayDifferAC-length(sigMicA)*L)/fs*v;
end

% if(DelayDifferAB<datanBC)
%     distDiffAB=(DelayDifferAB-1)/fs*v;
% else distDiffAB=(DelayDifferAB-length(sigMicA)*L-1)/fs*v;
% end

A=distDiffBC;
B=distDiffAC;
% C=distDiffAB;

%%判断是否在0度和180度，此时计算公式不能用%%%
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
    %%%计算A、B麦克传播距离r3/r2%%%
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
