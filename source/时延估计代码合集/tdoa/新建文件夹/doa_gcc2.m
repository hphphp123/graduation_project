function phi = doa_gcc2(alg,x1,x2,dx,N,Fs)
%function phi = doa_gcc2(alg,x1,x2,dx,N,Fs)
%
% direction estimation (azimuth phi) for 1 dim. microphone arrays
% using generalized cross correlation and speech activity detection
%
% (for use in oversampled filterbank)
%
% alg     1 SCOT-GCC
%         2 PHAT-GCC
% x1,x2   microphone signals
% dx      microphone distance in meters
% N       signal frame length (512, if omitted)
% Fs      sampling frquency in Hz (16000, if omitted)
% phi     azimuth in degrees
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

% algorithm: C.H.Knapp, G.C.Carter, "the generalized correlation method for 
%            estimation of time delay",
%            IEEE Trans ASSP, vol. ASSP-24, Aug. 1976, pp. 320-327. 

if nargin < 3
   help doa_gcc2
   return
end
if nargin < 5
   N = 512;
end
if nargin < 6
   Fs = 16000;
end

Lov = 4;
M = round(N/Lov);                            % frame hop size

% create time window function

a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1)));

x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2));
Ndata = 1+ceil((Nx-N)/M);
Nf = N;
Nfh = Nf/2+1;
Sxy = zeros(Nfh,1);
Sx = zeros(Nfh,1);
Sy = zeros(Nfh,1);
delay = zeros(Ndata,1);

dx = abs(dx);
OV = 4;                                      % oversampling factor for interpolation filter
vs = 340;                                    % acoustic waves propagation speed
Nd = 2+ceil(dx/vs*Fs);                       % max. delay (delay offset to obtain overall positive
                                             % delays in Cxy)
Nfo = OV*Nf;
Ndo = OV*Nd;

L = 2*Nd;

alpha = 0.8;                                % forgetting factor of spectral power averaging
alpha1 = 1-alpha;
doa_threshold = 0.01;                        % speech activity threshold
                                             % CHANGE, if necessary
delay_old = Ndo;

Cmat = zeros(L*OV,Ndata);

m = 0;
if alg == 1                                  % SCOT-GCC algorithm
  for n = 1:M:Nx-N+1
    m = m+1;
    n1 = n:n+N-1; 
    X1 = fft(x1(n1).*w,Nf);
    X2 = fft(x2(n1).*w,Nf);
    X1 = X1(1:Nfh);                          % use half of the spectra (real-valued signals)
    X2 = X2(1:Nfh);
    Sxy = alpha*Sxy + alpha1*X1 .* conj(X2); % spectra averaging
    Sx = alpha*Sx + alpha1*abs(X1).^2;
    Sy =  alpha*Sy + alpha1*abs(X2).^2;
    Cxy = OV*real(ifft(Sxy./(sqrt(Sx.*Sy)+1e-7),Nfo));  % generalized cross-correlation
    Cxy = [Cxy(Nfo-Ndo+1:Nfo) ; Cxy(1:Ndo)];
    Cmat(:,m) = Cxy; 
    [Cxymax,imax] = max(Cxy);
    if abs(Cxymax) > doa_threshold
       delay(m) = imax-1;                    % delay = maximum location of GCC
       delay_old = delay(m);
    else
       delay(m) = imax-1;
   end
  end
else                                         % PHAT-GCC algorithm
  for n = 1:M:Nx-N+1
    m = m+1;
    n1 = n:n+N-1; 
    X1 = fft(x1(n1).*w,Nf);
    X2 = fft(x2(n1).*w,Nf);
    X1 = X1(1:Nfh);
    X2 = X2(1:Nfh);
    Sxy = alpha*Sxy + alpha1*X1 .* conj(X2);
    Cxy = OV*real(ifft(Sxy./(abs(Sxy)+1e-4),Nfo));
    Cxy = [Cxy(Nfo-Ndo+1:Nfo) ; Cxy(1:Ndo)];
    Cmat(:,m) = Cxy; 
    [Cxymax,imax] = max(Cxy);
    if abs(Cxymax) > doa_threshold
       delay(m) = imax-1;
       delay_old = delay(m);
    else
       delay(m) = delay_old;
    end
  end
end
delay = mode(delay);
delay = (delay-Ndo)/(OV*Fs)            % correct delay offset
phi = 180/pi*real(acos(vs/dx*delay))
% 
% close all
% figurebackcolor = 'black';
% pos = [0.01 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
% 	     'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% t = M/Fs*[0:length(delay)-1];
% plot(t,phi)
% xlabel('Time t in sec.'), ylabel('Azimuth \phi in deg.');
% title(['d = ', num2str(100*dx),' cm,  ', 'F_s = ', num2str(Fs), ' Hz']);
% grid on, 
% axis tight
% %axis([0 50 0 180]);
% %set(gca,'PlotBoxAspectRatio',[1 0.5 0.5]);
% 
% % plot Cxy
% 
% pos = [0.5055 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','Generalized cross correlation',...
% 	     'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% 
% tau = 1000*linspace(-Nd/Fs,Nd/Fs,L*OV);
% imagesc(t,tau,Cmat); colormap(cool), colorbar
% %map = colormap('gray');
% %colormap(flipud(map));
% %imagesc(t,tau,Cmat); % colorbar;
% set(gca,'YDir','normal');
% ylabel('Delay \tau in msec.');
% xlabel('Time t in sec.');
% if alg == 1
%   title('Generalized SCOT cross correlation function R_{x_1x_2}(\tau,t)');
% else
%   title('Generalized PHAT cross correlation function R_{x_1x_2}(\tau,t)');
% end
% %set(gca,'PlotBoxAspectRatio',[1 0.5 0.5]);



