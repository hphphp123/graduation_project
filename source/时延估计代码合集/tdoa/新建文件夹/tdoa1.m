load original.mat
x1=data(:,3);
x2=data(:,2);

N=1024;
Lov = 4;
M = round(N/Lov);    
d=60;
v=34000;
%加窗
a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1)));

%%%%%%%%%%%

Nx=length(x1);
Ndata=1+ceil((Nx-N)/M); %帧数
Nf=N; %傅里叶变换长度
Nfh=Nf/2+1;  %频谱共轭对称，取一半
delay=zeros(Ndata,1); %每帧延迟矩阵
alpha = 0.8;                                % forgetting factor of spectral power averaging
alpha1 = 1-alpha;
d=abs(d); %距离为d
OV = 4; %内插倍数
Nd = 12;%最大延迟点数

Nfo = OV*Nf;%加倍的傅里叶变换长度
Ndo = OV*Nd;%加倍的最大延迟距离
L = 2*Nd;
alpha = 0.8;                                % forgetting factor of spectral power averaging
alpha1 = 1-alpha;
Sxy = zeros(Nf,1);
Sx = zeros(Nf,1);
Sy = zeros(Nf,1);
Cmat = zeros(L*OV,Ndata);

m = 0;
for n = 1:M:Nx-N+1
  m = m+1;
  n1 = n:n+N-1; 
  X1 = fft(x1(n1).*w,Nf);
  X2 = fft(x2(n1).*w,Nf);
  X1 = X1(1:Nf);                          
  X2 = X2(1:Nf);
  Sxy =alpha*Sxy + alpha1* X1 .* conj(X2); % spectra averaging
  Sx =alpha*Sx + alpha1* X1 .* conj(X1);
  Sy =  alpha*Sy + alpha1*X2 .* conj(X2);
  Cxy = OV*ifft(Sxy./(sqrt(abs(Sx.*Sy))),Nf);  % generalized cross-correlation 
%   Cxy = [Cxy(Nf-Nd+1:Nf) ; Cxy(1:Nd)];   %只能在前和后找，要不然就不对了
  [Cxymax,imax] = max(Cxy);
%   if imax<Nd
%     delay(m) = Nd-imax+1;
%   else 
      delay(m)=imax-1;
%   end
end
distDiffAC=mode(delay);