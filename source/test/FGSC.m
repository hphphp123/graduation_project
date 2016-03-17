% 
clear;
clc

M=4;   % 麦克风数目

clean = wavread('clean.wav');
x1 = wavread('noisy0.wav');
x2 = wavread('noisy10.wav');
x3 = wavread('noisy15.wav');
x4 = wavread('noisy20.wav');
% x5 = wavread('noisy25.wav');
Len = length(clean);
mic1 = x1(1:Len);
mic2 = x2(1:Len);
mic3 = x3(1:Len);
mic4 = x4(1:Len);
% mic5 = x5(1:Len);
x = [mic1';mic2';mic3';mic4'];%;mic5'
FBFout=sum(x)/M;

B = [1 -1 0 0 
     0 1 -1 0 
     0 0 1 -1 ];
Bout = B*x;
x1 = Bout(1,:)';
x2 = Bout(2,:)';
x3 = Bout(3,:)';
% x4 = Bout(4,:)';

y1 = zeros(1,Len);
y2 = y1;
y3 = y1;
y4 = y1;
MCout = y1;
N = 64;
h1 = zeros(1,N);
h2 = h1;
h3 = h1; 
h4 = h1;
GSCout = zeros(1,Len);
u = 0.0002;
weight = [];
for i=1:fix(Len/N)-1
    X1 = fft(x1((i-1)*N+1:(i+1)*N));
    H1 = fft([h1,zeros(1,N)]);
    O11 = real(ifft(X1'.*H1));
    y1(i*N+1:(i+1)*N) = O11(N+1:N*2);
    
    X2 = fft(x2((i-1)*N+1:(i+1)*N));
    H2 = fft([h2,zeros(1,N)]);
    O12 = real(ifft(X2'.*H2));
    y2(i*N+1:(i+1)*N) = O12(N+1:N*2);
    
    X3 = fft(x3((i-1)*N+1:(i+1)*N));
    H3 = fft([h3,zeros(1,N)]);
    O13 = real(ifft(X3'.*H3));    
    y3(i*N+1:(i+1)*N) = O13(N+1:N*2);
    
%     X4 = fft(x4((i-1)*N+1:(i+1)*N));
%     H4 = fft([h4,zeros(1,N)]);
%     O14 = real(ifft(X4'.*H4));
%     y4(i*N+1:(i+1)*N) = O14(N+1:N*2);
    
    X = [X1';X2';X3'];%;X4'
    MCin = sum(X); 
    MCout(i*N+1:(i+1)*N) = sum([y1(i*N+1:(i+1)*N);y2(i*N+1:(i+1)*N);y3(i*N+1:(i+1)*N)]);%;y4(i*N+1:(i+1)*N)
    
    GSCout(i*N+1:(i+1)*N) = FBFout(i*N+1:(i+1)*N)-(MCout(i*N+1:(i+1)*N));
    E = fft([zeros(1,N),GSCout(i*N+1:(i+1)*N)]);
    O2 = real(ifft(E.*conj(MCin)));
    V = O2(1:N);
    
    h1 = h1+2*u*V;
    h2 = h1+2*u*V;
    h3 = h3+2*u*V;
    %h4 = h4+2*u*V;
    
    weight = [weight h1(1)];
end
subplot(211);plot(mic1);wavplay(mic1)
subplot(212);plot(GSCout);wavplay(GSCout)
snr2=SNR_singlech(clean,mic1);fprintf('snr2=%5.1f\n',snr2);
snr1=SNR_singlech(clean,GSCout);fprintf('snr1=%5.1f\n',snr1);
%wavwrite(GSCout,'Fgsc.wav');