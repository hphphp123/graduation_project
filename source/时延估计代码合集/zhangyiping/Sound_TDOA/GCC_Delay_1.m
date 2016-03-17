
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



%%% ��֡+�ڲ壬��ʱ�ӹ���delay������
x1=sigMicA ;
x2=sigMicB;
x3=sigMicC;

N=2048;
Lov = 4;
M = round(N/Lov);    
%�Ӵ�
a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1)));

Nx=length(x1);
Ndata=1+ceil((Nx-N)/M); %֡��
Nf=N;  %����Ҷ�任����
Nfh=Nf/2+1;   %Ƶ�׹���Գƣ�ȡһ��
Sxy=zeros(Nf,1);  %��������
Sx=zeros(Nf,1);  %������
Sy=zeros(Nf,1);
delay=zeros(Ndata,1);  %ÿ֡�ӳپ���
d=abs(d); %����Ϊd
Nd = 2+ceil(d/v*fs);%����ӳٵ���

alpha = 0.8;    %�������ʣ������źŷ�ƽ��                          
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
  Cxy = [Cxy(Nf-Nd+1:Nf) ; Cxy(1:Nd)];   %ֻ����ǰ�ͺ��ң�Ҫ��Ȼ�Ͳ�����
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
  Cxy = [Cxy(Nf-Nd+1:Nf) ; Cxy(1:Nd)];   %ֻ����ǰ�ͺ��ң�Ҫ��Ȼ�Ͳ�����
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
