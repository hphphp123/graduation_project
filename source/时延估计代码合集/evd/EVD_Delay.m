%% lms
function [thta,A,B] = EVD_Delay(data ,parameter)

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



Nd = 2+ceil(d/v*fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Ndo = Nd*Lov;
mu=0.3;

doa_threshold = -0.09;    % speech activity threshold
                          % CHANGE, if necessary

Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)
N=512;
M=3500;
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
Sx3 = zeros(Nf,1);
alpha = 0.2;              % forgetting factor of spectral power averaging  
alpha1 = 1-alpha;

% loop to compute eigenvector u corresponging to zero eigenvalue

k = 0;
mb = 0;
for n = 1:N:Nx-Nf
   mb = mb+1;
   m = n:n+Nf-1;
   X1 = fft(x1(m));
   X3 = fft(x3(m));
   e = real(ifft(U1.*X1+U2.*X3));
   E = fft([zeros(N,1) ; e(N+1:Nf)]);
   Sx1 = alpha*Sx1 + alpha1*abs(X1).^2;
   Sx3 = alpha*Sx3 + alpha1*abs(X3).^2;
   U1 = U1 - (mu ./ (Sx1+eps)) .* conj(X1) .* E;
   U2 = U2 - (mu ./ (Sx3+eps)) .* conj(X3) .* E;
   if mod(mb,Mb) == 0           % find delay, and restart adaptive filter
      u2 = real(ifft(U2));      % eigenvector to be used to find delay
      u2 = u2(Nh-Nd:Nh+Nd-1);
      u2 = resample(u2,Lov,1);  % interpolate
      k = k+1;
      U(:,k) = u2;
      [umin,dmin] = min(u2);
      del = dmin-Ndo;
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
delay_evd1=delay(end)/(fs*Lov)*v;
Nd = 2+ceil(d/v*fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Ndo = Nd*Lov;
mu=0.3;

doa_threshold = -0.09;    % speech activity threshold
                          % CHANGE, if necessary

Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)
N=512;
M=3500;
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
Sx3 = zeros(Nf,1);
alpha = 0.2;              % forgetting factor of spectral power averaging  
alpha1 = 1-alpha;

k = 0;
mb = 0;
for n = 1:N:Nx-Nf
   mb = mb+1;
   m = n:n+Nf-1;
   X2 = fft(x2(m));
   X3 = fft(x3(m));
   e = real(ifft(U1.*X2+U2.*X3));
   E = fft([zeros(N,1) ; e(N+1:Nf)]);
   Sx2 = alpha*Sx2 + alpha1*abs(X2).^2;
   Sx3 = alpha*Sx3 + alpha1*abs(X3).^2;
   U1 = U1 - (mu ./ (Sx2+eps)) .* conj(X2) .* E;
   U2 = U2 - (mu ./ (Sx3+eps)) .* conj(X3) .* E;
   if mod(mb,Mb) == 0           % find delay, and restart adaptive filter
      u2 = real(ifft(U2));      % eigenvector to be used to find delay
      u2 = u2(Nh-Nd:Nh+Nd-1);
      u2 = resample(u2,Lov,1);  % interpolate
      k = k+1;
      U(:,k) = u2;
      [umin,dmin] = min(u2);
      del = dmin-Ndo;
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
delay_evd2=delay(end)/(fs*Lov)*v;


A=delay_evd2;
B=delay_evd1;
M=2*A^2+B^2-4*A*B+2*d^2;
N=-2*(2*A-B);
r3=M/N;
r2=A+r3-B;

cosA=(r3^2+d^2-r2^2)/(2*d*r3);
cosA = min(cosA,1);
thta = acos(cosA)*180/pi;


