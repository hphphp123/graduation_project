clc;clear;
%% ������FFT�������ֱ��ʹ�ã������н�ȡ 
% [ns,fs,nbits]=wavread('E:\speech enhancement\ANR_Process\ANR_Process\test\sp10_airport_sn0.wav');
fid=fopen('E:\������ǿMATLAB\NOIZEUS\shiyanshinoise3.pcm','r');
ns00=fread(fid,'int16');
ns=[zeros(160,1);ns00];
fs=8000;

% x=ns00(1:320,1);
% [y1,y2,pfY1,pfY2]=smoothing (x,1);
% Y1=pfY1';

ninc=round(0.020*fs); %%��Ӧ����Ƶ��֡��
if rem(ninc,2)==1, ninc=ninc+1; end;
ovf=2;                  % �ص����ӣ�50%����
Len=ninc*ovf;           % ���η�������   
Len1=Len-ninc;
%win=hamming(Len,'periodic');%���崰
for swCnt = 0:1:(Len-1)
	win((swCnt+1),1) = 0.5-0.5*cos(2*pi*swCnt/Len);
end
    
k=1;  %��ǰ����֡��ʼ�����
% nFFT=2^nextpow2(Len);
nFFT=Len;
Nframes=floor(length(ns)/ninc)-1;
img=sqrt(-1);
x_old=zeros(Len1,1);
xfinal=zeros(Nframes*ninc,1);

for n=1:Nframes
    ns_cut=ns(k:k+Len-1); 
    ns_win=win.*ns(k:k+Len-1);
    ns_fft=fft(ns_win,nFFT)/nFFT;
%     ns_out=ns_fft;
    ns_ps1=ns_fft.*conj(ns_fft); %������������
    
    line(:,1)= ns_ps1(1:160);
    ns_ps(:,1)=VOICE_FilterLineToBand(line);
%     ns_cutN(:,n)=ns_cut; 
%     ns_WinN(:,n)=ns_win; 
%     ns_fftN(:,n)=ns_fft; 
%     ns_psN(:,n)=ns_ps1; 
      ns_psB(:,n)=ns_ps; 
    
    if n == 1
      [parametersN,parametersE] = imcra_initialise_parameters(ns_ps); %���ݳ�ʼ�� 
    elseif n<=60
      
      ns_ps1 = floor(ns_ps);
      parametersN = estnoise_imcra(ns_ps1,parametersN);
      parametersE.noise_ps=parametersN.noise_ps;
      aa=parametersE.noise_ps*1024;
      parametersE = logmmse_OM(ns_ps,parametersE);
      parametersN.GH1=parametersE.GH1;
      
    else
      parametersN = estnoise_imcra(ns_ps,parametersN);
      parametersE.noise_ps=parametersN.noise_ps;
      parametersE = logmmse_OM(ns_ps,parametersE);
      parametersN.GH1=parametersE.GH1;  
    end
    
    noise_ps(:,n)=parametersN.noise_ps; 
    SminB(:,n)=parametersN.Smin;  
    SfB(:,n)=parametersN.Sf; 
    SB(:,n)=parametersN.S; 
    GH1B(:,n)=parametersN.GH1; 
    gammaB(:,n)=parametersN.gamma*512; 

    
    sig1 =  parametersE.sigout; 
    line=VOICE_FilterBandToLin(sig1);
    sig(1:160,1)=line(1:160);
    sig(161,1)=0;
    sig(162:320,1)=flipud(line(2:160,1));

    xi_w= ifft(sig .* (ns_fft*nFFT));
    xi_w= real( xi_w);

% --------- Overlap and add ---------------
%
    xfinal(k:k+ninc-1)= x_old+ xi_w(1:Len1);
    x_old= xi_w(Len1+1:Len);  
    k=k+ninc;
end

ns16=int16(xfinal);
fid=fopen('E:\������ǿMATLAB\shiyanshinoise3_out.pcm','w');
fwrite(fid,ns16,'int16');
%wavwrite(xfinal,fs,16,'E:\speech
%enhancement\ANR_Process\ANR_Process\test\capture_speech_KFC_out1.wav');
