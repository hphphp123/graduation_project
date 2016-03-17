%基于麦克风小阵的多噪声环境语音增强算法，
%关键词：非平稳噪声：麦克风阵列：广义旁瓣抵消：相干滤波语音增强

close all;clc;clear all;
tic;
[s,fs,bits]=wavread('d:\语音文件\clean\sp01.wav');                 %纯语音信号
s=s-mean(s);                   %消除直流分量
s=s/max(abs(s));               %幅值归一化
N=length(s);                   %语音信号长度
noise1=wavread('d:\noisex-92\white.wav');
[x1,noise]=add_noisedata(s,noise1,fs,fs,5); %第一路添加白噪的带噪语音
k=5;
x1=x1';
%r=zeros(1,2000);
r=wgn(1,N,-20);%产生-20dB高斯噪声
% r=randn(1,2000); %noise
% r=r/max(abs(r));
x2=[r(1:k) x1(1:N-k)] ;          %delay signal
noise2=wavread('d:\noisex-92\babble.wav');
[x2,noise]=add_noisedata(x2,noise2,fs,fs,5); %第二路添加白噪的带噪语音
n=1:N;
% figure 
% subplot(211)
% plot(n,x1);
% axis([1 N -1 1]);
% title('模拟第一路麦克风信号');
% subplot(212)
% plot(n,x2);
% axis([1 N -1 1]);
% title('延迟后模拟第二路麦克风信号')
%% fastlms

vs = 340;                 % acoustic waves propagation speed
dx=0.1;
dx = abs(dx);
Nd = 2+ceil(dx/vs*fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Fs1 = Lov*fs;
Ndo = Nd*Lov;
mu=0.3;
N=512;
M=2048;

doa_threshold = 0.2;      % speech activity threshold
                          % CHANGE, if necessary

x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)

Mb = max(1,round(M/N));
M = Mb*N;
Nf = 2*N;                 % FFT length

if Nx <= Nf
   disp('INFO: FFT length must be less than signal length');
   return
end

% init. variables

Nh = floor(N/2);          % delay of first microphone signal
W = zeros(Nf,1);          % weight vector in frequency domain
Nb = ceil((Nx-Nf)/N);
delay = zeros(Nb,1);
delay_old = 0;
Wmat = zeros(2*Ndo,Nb);

Sx = zeros(Nf,1);         % spectral power used with step size mu
alpha = 0.2;              % forgetting factor of spectral power averaging
alpha1 = 1-alpha;

% adaptive filter loop

k = 0;
mb = 0;
for n = 1:N:Nx-Nf
   mb = mb+1;
   X2 = fft(x2(n:n+Nf-1));
   y = real(ifft(W.*X2));
   e = x1(n+Nh+1:n+Nh+N) - y(N+1:Nf);
   E = fft([zeros(N,1) ; e]);
   Sx = alpha*Sx + alpha1*abs(X2).^2;
   W = W + (mu ./ (Sx + eps)) .* conj(X2) .* E; 
   if mod(mb,Mb) == 0          % determine delay every Mb blocks only
      w = real(ifft(W));
      w1 = w(Nh-Nd:Nh+Nd-1);   % weight vector to be used to find delay
      w1 = resample(w1,Lov,1); % interpolate
      k = k+1;
      Wmat(:,k) = w1;          % filter coefficient map
      [wmax,dmax] = max(w1);
      del = abs(dmax-Ndo)+1;
      if wmax < doa_threshold  
         delay(k) = delay_old; % use old delay during speech pauses 
      else
         delay(k) = del;
         delay_old = del;
      end
   end
end


delay = delay(1:k);
disp(sprintf('final delay = %2.4e msec (%d samples)', ...
             1000*delay(end)/Fs1, delay(end)));
phi = 180/pi*real(acos(vs/dx*delay/Fs1));
disp(sprintf('final azimuth = %2.4e deg', phi(end)));
%% 子空间分解
vs = 340;                 % acoustic waves propagation speed
dx=0.1;
dx = abs(dx);
Nd = 2+ceil(dx/vs*fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Fs1 = Lov*fs;
Ndo = Nd*Lov;
N=512;
M=3500;
mu=0.3;

doa_threshold = -0.09;    % speech activity threshold
                          % CHANGE, if necessary

x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)

Mb = max(1,round(M/N));
M = Mb*N;
Nf = 2*N;                 % FFT length   

% init. eigenvector u

Nh = floor(N/2)+1;
u0_1 = zeros(N,1);
u0_1(Nh) = 1;
U0_1 = fft(u0_1,Nf);
U0_2 = zeros(Nf,1);
U1 = U0_1;
U2 = U0_2;
Nb = ceil((Nx-N)/M);
delay = zeros(Nb,1);
delay_old = 0;
U = zeros(2*Ndo,Nb);

Sx1 = zeros(Nf,1);        % spectral power used with step size mu
Sx2 = zeros(Nf,1);
alpha = 0.2;              % forgetting factor of spectral power averaging  
alpha1 = 1-alpha;

% loop to compute eigenvector u corresponging to zero eigenvalue

k = 0;
mb = 0;
for n = 1:N:Nx-Nf
   mb = mb+1;
   m = n:n+Nf-1;
   X1 = fft(x1(m));
   X2 = fft(x2(m));
   e = real(ifft(U1.*X1+U2.*X2));
   E = fft([zeros(N,1) ; e(N+1:Nf)]);
   Sx1 = alpha*Sx1 + alpha1*abs(X1).^2;
   Sx2 = alpha*Sx2 + alpha1*abs(X2).^2;
   U1 = U1 - (mu ./ (Sx1+eps)) .* conj(X1) .* E;
   U2 = U2 - (mu ./ (Sx2+eps)) .* conj(X2) .* E;
   if mod(mb,Mb) == 0           % find delay, and restart adaptive filter
      u2 = real(ifft(U2));      % eigenvector to be used to find delay
      u2 = u2(Nh-Nd:Nh+Nd-1);
      u2 = resample(u2,Lov,1);  % interpolate
      k = k+1;
      U(:,k) = u2;
      [umin,dmin] = min(u2);
      del = abs(dmin-Ndo)+1;
      if umin > doa_threshold   % signal to weak (e.g. speech pauses)
         delay(k) = delay_old;
      else
         delay(k) = del;
         delay_old = del;
      end
      U1 = U0_1;
      U2 = U0_2;
   end
end

delay = delay(1:k);
disp(sprintf('final delay = %2.4e msec (%d samples)', ...
             1000*delay(end)/Fs1, delay(end)));
phi = 180/pi*real(acos(vs/dx*delay/Fs1));
disp(sprintf('final azimuth = %2.4e deg', phi(end)));
%%
y=xcorr(x2,x1);    %Cross-correlation of the signal from-130 to 130
x=find(y==max(y)) ;            %延迟
delay_cc=x-length(x1)
% subplot(313)
% x2=x2';
% x3=[x2(delay:N) r(1:delay-1)];
% plot(n,x3);
% axis([1 N -1 1]);
% title('时延估计对齐后得到的第二路麦克风信号');
x1=x1';
% x2=x2';
s=s';
% snr1=SNR_singlech(s(1:N),x1);fprintf(' 第一路麦克风信号snr1=%5.1f\n',snr1); 
% snr2=SNR_singlech(s(1:N),x2);fprintf(' 第二路麦克风信号snr2=%5.1f\n',snr2);
% sound(x1);
% pause(1);
% sound(x2);
X1=fft(x1,length(x1));
X2=fft(x2',length(x2));
FXX=X1.*conj(X2);
Amplitude=abs(X1).*abs(X2);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_phat=length(x1)-DelayDifferAB+1
X1=fft(x1,length(x1));
X2=fft(x2,length(x2));
FXX1=X1.*conj(X1);
FXX2=X2.*conj(X2);
Amplitude=sqrt(abs(FXX1.*FXX2));
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_scot=length(x1)-DelayDifferAB+1
X1=fft(x1,length(x1));
X2=fft(x2,length(x2));
FXX1=X1.*conj(X1);
FXX2=X2.*conj(X2);
Amplitude=abs(FXX1);
GCC=FXX./Amplitude;
R_XX=ifft(GCC);
[val1,DelayDifferAB]=max(R_XX);
delay_roth=length(x1)-DelayDifferAB+1
len=length(x1);
x=[zeros(1000,1);x1];
w=zeros(1000,1);%初始1000阶加权系数
u=0.000005;
y=zeros(len,1);
for i=1:len-4000
     y(i)=(x(i+999:-1:i))'*w;     
     e=x2(i)-y(i); 
     w=w+2*u*e*x(i+999:-1:i);
end
% subplot(411),plot(x);grid on
% subplot(412),plot(x2);grid on
% subplot(413),plot(y);grid on
% subplot(414),plot(w);grid on
%wavwrite(y,'处理后信号');
y1=xcorr(x,x2);
%subplot(515),plot(y1);grid on
%delay_lms=find(y1==max(y1))-len
delay_lms=find(w==max(w))