function y = adaptive_beam(Xin,mics,phi_look,mu,Gw,f,Fs)
%function y = adaptive_beam(X,mics,phi_look,mu,Gw,f,Fs)
%
% adaptive beamformer using an FFT filter bank with/without DOA estimation and tracking
% (frequency domain version of Frost beamformer)
% version to plot array pattern (image and 3d plot)
%
% X          matrix of microphone signals (one column for each channel)
% mics       microphone placement
%            (first col = x, second col = y coord.)
% phi_look   azimuth in deg. of array look direction (0 ... 180°)
%            if phi_look=[] (or omitted), then automatic speaker tracking is used 
% mu         step size of adaptive algorithm (default 0.01)
% Gw         Gw parameter in dB limit at low frequencies (default 4)
% f          frequency of array pattern in Hz (default 1000 Hz) 
% Fs         sampling frequency in Hz (default 16000 Hz)
% y          beamformer output signal
%
% example:   y = adaptive_beam(X,mics);    % automatic speaker tracking
%            y = adaptive_beam(X,mics,90); % look direction = 90°    
%
%   Copyright 2006 Gerhard Doblinger, Vienna University of Technology
%   g.doblinger@tuwien.ac.at
%   http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


if nargin == 0
   help adaptive_beam
   return
elseif nargin == 2
   phi_look = [];
   mu = 0.01;
   Gw = 4;
   f = 1000;
   Fs = 16000;
elseif nargin == 3
   mu = 0.01;
   Gw = 4;
   f = 1000;
   Fs = 16000;
elseif nargin == 4
   Gw = 4;
   f = 1000;
   Fs = 16000;
elseif nargin == 5
   f = 1000;
   Fs = 16000;
elseif nargin == 6
   Fs = 16000;
end

if f >= Fs/2
   error('f >= Fs/2');
end

[Nx,N] = size(Xin);
[Nc,dum] = size(mics);
if Nc ~= N
   error('number of mics does not match number of input signals');
end

Nfft = 512;                          % FFT length = time window length
Nfh = Nfft/2+1;
M = Nfft/4;                          % frame hop size 
disp(sprintf('Nmics = %d, Nfft = %d, Fs = %5.0f Hz',N, Nfft, Fs));

% set max. input signal frequency according to min. array spacing

vs = 340;                            % acoustic waves propagation speed
dmin = min(diff(mics(:,1)));
fmax = min(Fs/2,vs/(2*dmin));
f = min(fmax,f);
Nfmax = round(Nfft*fmax/Fs);
disp(sprintf('Nfmax = %d, fmax = %5.0f Hz', Nfmax, fmax));

% compute parameters of starting solution (at azimuth phi_d)  

theta_d = 90;
if ~isempty(phi_look)
   phi_look = abs(phi_look);
   phi_d = phi_look;
else
   phi_d = 90;
end
[b,Wc,P] = compute_params(mics,vs,Nfft,Nfmax,theta_d,phi_d,Gw,Fs);

% variables of DOA estimation (GCC algorithms)

mic1 = 1;                            % index of first and second mic
mic2 = N;
doa_threshold = 0.22;                % CHANGE, if necessary!
Sxy = zeros(1,Nfh);
Ndata = 1+ceil((Nx-Nfft)/M);
if ~isempty(phi_look)
   phi_doa = phi_look*ones(Ndata,1);
   disp(sprintf('no speaker tracking, array look direction = %d deg.', phi_look));
else
   disp('automatic speaker tracking enabled');
   phi_doa = zeros(Ndata,1);
end

OV = 4;                              % oversampling factor for Cxy interpolation
dx = abs(mics(mic1,1)-mics(mic2,1)); % distance of microphones
Nd = 2+ceil(dx/vs*Fs);               % max. delay (delay offset to obtain overall positive
                                     % delays in Cxy)
Nfov = OV*Nfft;
Ndov = OV*Nd;
Tsov = 1/(OV*Fs);
delay_old = 0;                       % starting delay, and azimuth
phi_old = 90;                        % (used if tracking enabled)

alpha = 0.8;                         % forgetting factor of Sxy averaging
alpha1 = 1-alpha;

% variabels of array pattern

phi = pi/180*[0:180];
Nphi = length(phi);
er = mics*sin(pi/180*theta_d(1))*[cos(phi) ; sin(phi)];
mf1 = round(Nfft*f/Fs);              % nearest frequency index
beta = 2*pi*mf1*Fs/(vs*Nfft);
D1 = exp(j*beta*er);                 % matrix of steering vectors
Ndata = 1+ceil((Nx-Nfft)/M);
RdB = zeros(Ndata,Nphi);

% filter bank 

W = Wc;                              % initial adaptive filter coefficients
V = zeros(size(Wc));
y = zeros(Nx,1);
h = hanning(Nfft);
H = h(:) * ones(1,N);                % matrix of time windows

disp('computing beam pattern, and adaptive array output...');
drawnow

n = 0;
for  m = 1:M:Nx-Nfft+1
   
   % compute filter bank output signal
   
   n = n+1;
   m1 = m:m+Nfft-1; 
   X = fft(Xin(m1,:) .* H).';        % apply FFT to weighted frames of all mics
   Y = sum(conj(W) .* X(:,1:Nfmax)); % apply adaptive filter coefficients, sum up FFT bins
   yb = real(ifft(Y,Nfft));
   y(m1) = y(m1) + yb(:);            % overlap-add
   
   if isempty(phi_look)              % DOA estimation (PHAT-GCC)
      X1 = X(mic1,1:Nfh);
      X2 = X(mic2,1:Nfh);
      Sxy = alpha*Sxy + alpha1*X1 .* conj(X2);     % spectra averaging
      Cxy = OV*real(ifft(Sxy./(abs(Sxy)+1e-4),Nfov));
      Cxy = [Cxy(Nfov-Ndov+1:Nfov) Cxy(1:Ndov)];
      [Cxymax,imax] = max(Cxy);
      if Cxymax > doa_threshold
         delay = (imax-1-Ndov)*Tsov; % delay = maximum location of GCC
         delay_old = delay;
         factiv = 1;
      else
         delay = delay_old;
         factiv = 0;
      end
      phi_doa(n) = 180/pi*real(acos(vs/dx*delay));   
   
      % compute new starting solution for new estimated direction
   
      if factiv                      % speech activity?
         if abs(phi_doa(n)-phi_old) >= 2           % new start with new parameters  
            [b,Wc,P] = compute_params(mics,vs,Nfft,Nfmax,theta_d,phi_doa(n),Gw,Fs);
            phi_old = phi_doa(n);    % save DOA estimate
         end
      end
   end
   
   % update coefficient matrix W (constrained LMS algorithm)
   
   for k = 1:Nfmax
     V(:,k) = P(:,:,k) * (V(:,k) - mu/(norm(X(:,k))^2+eps)*X(:,k)*conj(Y(k)));
     no = norm(V(:,k));
     if no > b(k)                    % norm constraint
       V(:,k) = V(:,k)*b(k)/no;
     end
   end
   W = Wc + V;
   
   % compute array pattern

   R = abs(W(:,mf1)'*D1).^2;
   RdB(n,:) = R/max(R);

end

y = y/Nfft*M*4;

phi_doa = phi_doa(1:n);

% plot 2-dim. array pattern

close all
pos = [0.01 0.5 0.49 0.4];
figurebackcolor = 'black';
fabf = figure('numbertitle','off','name','Adaptive beamformer I',...
	     'Units','normal','Position',pos, 'DoubleBuffer', 'on');
colordef(fabf,figurebackcolor);
RdB = RdB(1:n,:);
RdB = RdB';
RdB = max(-40,10*log10(RdB+eps));
t = M/Fs*[0:n-1];
map = colormap('jet');
%colormap(flipud(map));

hold on
grid on
imagesc(t,180/pi*phi,RdB);
colorbar
plot(t,phi_doa, 'y');
hold off
set(gca,'YDir','normal');
ylabel('Azimuth \Phi in deg.');
xlabel('Time t in sec.');
title(['Array pattern in dB at f = ', num2str(f),' Hz, \mu = ', num2str(mu)]);
text(1.04*t(end),185,'  dB');
if isempty(phi_look)
   text(3,96,'Estimated Azimuth', 'Color', 'yellow');
else
   text(3,phi_look+6,'Look Direction', 'Color', 'yellow');
end
axis([0,t(end),0,180]);

% plot 3-dim. array pattern

pos(1) = pos(1)+0.5;
   fp = figure('numbertitle','off','name','Adaptive beamformer II',...
	       'Units','normal','Position',pos);
colordef(fp,figurebackcolor);

id = 1:50:length(t);
colormap('gray');
surfc(t,180/pi*phi,RdB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
surfl(t(id),180/pi*phi,RdB(:,id));
shading interp
ylabel('Azimuth \Phi in deg.');
xlabel('Time t in sec.');
zlabel('dB');
title(['Array pattern in dB at f = ', num2str(f),' Hz, \mu = ', num2str(mu)]);
axis([0,t(end),0,180,-40,0]);
view([30,45]);


%---------------------------------------------------------------------------

function w = delaysum_beam(mics,vs,theta_d,phi_d,f);
% function w = delaysum_beam(mics,theta_d,phi_d,f);
%
% compute weights of delay&sum beamformer
%
% w              weight vector of length N (number of mics)
% mics           [xi,yi] microphone cooerdinates (N x 2 matrix)
% vs             speed of spound im m/sec  
% theta_d        desired elevation in deg.
% phi_d          desired azimuth in deg.
% f              frequency in Hz
%
% G.Doblinger, TU-Wien 10-04

if nargin ~= 5
   help delaysum_beam;
   return;
end

theta_d = theta_d*pi/180;
phi_d = phi_d*pi/180;

beta = 2*pi*f/vs;         % wave number
[N,K] = size(mics);
if (K < 2) | (N < 1)
   error('bad matrix of microphine coordinates');
end

if K == 2
   rn = [mics zeros(N,1)];
else
   rn = mics;
end

% steering vector of desired direction

ed = [sin(theta_d)*cos(phi_d); sin(theta_d)*sin(phi_d); cos(theta_d)];
d0 = exp(j * beta * rn * ed);
w = d0 / N;

%---------------------------------------------------------------------------

function [b,Wc,P] = compute_params(mics,vs,Nfft,Nfmax,theta_d,phi_d,Gw,Fs)

% compute parameters of starting solution

[N,dum] = size(mics);                % N = number of microphones
if Nfmax >  Nfft/2+1
   error('Nfmax  > Nfft/2+1');
end
Nco = length(theta_d);               % number of constraints
P = zeros(N,N,Nfmax);
Wc = zeros(N,Nfmax);
b = zeros(Nfmax,1);
if Nco == 1
   for l = 1:Nfmax
      if ((l-1)/Nfft*Fs > 250) & (Gw < -8)  Gw = -8; end
      if ((l-1)/Nfft*Fs > 450) & (Gw < -2)  Gw =  -2; end
      if ((l-1)/Nfft*Fs > 700) & (Gw < 2)  Gw =  2; end
      if ((l-1)/Nfft*Fs > 1000) & (Gw < 4) Gw =  4; end
      if ((l-1)/Nfft*Fs > 2000) & (Gw < 6) Gw =  6; end
      c = delaysum_beam(mics,vs,theta_d,phi_d,(l-1)/Nfft*Fs);
      P(:,:,l) = eye(N) - N*c*c';
      Wc(:,l) = c;
      b(l) = sqrt(10^(-Gw/10) - Wc(:,l)' * Wc(:,l));     
   end
else
   g = [1;zeros(Nco-1,1)];           % vector of constraints
   C = zeros(N,Nco,Nfmax);
   Creg = 1e-10*eye(Nco);            % add regularization term to avoid singular matrix
   for l = 1:Nfmax
      if ((l-1)/Nfft*Fs > 250) & (Gw < -8)  Gw = -8; end
      if ((l-1)/Nfft*Fs > 450) & (Gw < -2)  Gw =  -2; end
      if ((l-1)/Nfft*Fs > 700) & (Gw < 2)  Gw =  2; end
      if ((l-1)/Nfft*Fs > 1000) & (Gw < 4) Gw =  4; end
      if ((l-1)/Nfft*Fs > 2000) & (Gw < 6) Gw =  6; end
      for co = 1:Nco
	 C(:,co,l) = N*delaysum_beam(mics,vs,theta_d(co),phi_d(co),(l-1)/Nfft*Fs);
      end
      C1 = C(:,:,l) * (Creg + C(:,:,l)' * C(:,:,l))^-1; 
      P(:,:,l) = eye(N) - C1*C(:,:,l)';
      Wc(:,l) = C1 * g;
      b(l) = sqrt(10^(-Gw/10) - Wc(:,l)' * Wc(:,l));     
   end
end
b = real(b);                         % remove small imaginary part due to rounding

%---------------------------------------------------------------------------




