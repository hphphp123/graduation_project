%% lms
function [thta,A,B] = LMS_Delay(data ,parameter)

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

x1=sigMicA ;
x2=sigMicB;
x3=sigMicC;


d = abs(d);
Nd = 2+ceil(d/v*fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Fs1 = Lov*fs;
Ndo = Nd*Lov;

doa_threshold = 0.2;      % speech activity threshold
                          % CHANGE, if necessary

x1 = x1(:);
x2 = x2(:);
x3 = x3(:);

Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)
N=512;
M=4*N;
mu=0.25;
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
   X3 = fft(x3(n:n+Nf-1));
   y = real(ifft(W.*X3));
   e = x1(n+Nh+1:n+Nh+N) - y(N+1:Nf);
   E = fft([zeros(N,1) ; e]);
   Sx = alpha*Sx + alpha1*abs(X3).^2;
   W = W + (mu ./ (Sx + eps)) .* conj(X3) .* E; 
   if mod(mb,Mb) == 0          % determine delay every Mb blocks only
      w = real(ifft(W));
      w1 = w(Nh-Nd:Nh+Nd-1);   % weight vector to be used to find delay
      w1 = resample(w1,Lov,1); % interpolate
      k = k+1;
      Wmat(:,k) = w1;          % filter coefficient map
      [wmax,dmax] = max(w1);
      del = dmax-Ndo;
      if wmax < doa_threshold  
         delay(k) = delay_old; % use old delay during speech pauses 
      else
         delay(k) = del;
         delay_old = del;
      end
   end
end

delay = delay(1:k);
delay_lms1=delay(end)/(Lov*fs)*v;

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
   X3 = fft(x3(n:n+Nf-1));
   y = real(ifft(W.*X3));
   e = x2(n+Nh+1:n+Nh+N) - y(N+1:Nf);
   E = fft([zeros(N,1) ; e]);
   Sx = alpha*Sx + alpha1*abs(X3).^2;
   W = W + (mu ./ (Sx + eps)) .* conj(X3) .* E; 
   if mod(mb,Mb) == 0          % determine delay every Mb blocks only
      w = real(ifft(W));
      w1 = w(Nh-Nd:Nh+Nd-1);   % weight vector to be used to find delay
      w1 = resample(w1,Lov,1); % interpolate
      k = k+1;
      Wmat(:,k) = w1;          % filter coefficient map
      [wmax,dmax] = max(w1);
      del = dmax-Ndo;
      if wmax < doa_threshold  
         delay(k) = delay_old; % use old delay during speech pauses 
      else
         delay(k) = del;
         delay_old = del;
      end
   end
end

delay = delay(1:k);
delay_lms2=delay(end)/(Lov*fs)*v;


A=delay_lms2;
B=delay_lms1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta = acos(cosA)*180/pi;


