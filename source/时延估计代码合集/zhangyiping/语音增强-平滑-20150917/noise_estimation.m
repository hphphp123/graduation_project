clc;clear;
%% 程序中FFT分析结果直接使用，不进行截取 
[ns00,fs,nbits]=wavread('E:\语音增强MATLAB\NOIZEUS\00_ns_road_ori.wav');
% [nn,fs,nbits]=wavread('E:\语音增强MATLAB\noise.wav');
ns=[zeros(160,1);ns00];
ninc=round(0.020*fs); %%对应采样频率帧长
if rem(ninc,2)==1, ninc=ninc+1; end;
ovf=2;                  % 重叠因子，50%交叠
Len=ninc*ovf;           % 单次分析长度   
Len1=Len-ninc;
win=hamming(Len,'periodic');%定义窗
k=1;  %当前分析帧起始点序号
% nFFT=2^nextpow2(Len);
nFFT=Len;
Nframes=floor(length(ns)/ninc)-1;
img=sqrt(-1);
x_old=zeros(Len1,1);
xfinal=zeros(Nframes*ninc,1);

for n=1:Nframes
    ns_cut=ns(k:k+Len-1); 
    ns_win=win.*ns(k:k+Len-1);
    ns_fft=fft(ns_win,nFFT);
    ns_out(:,n)=abs(ns_fft);
    ns_ps=ns_fft.*conj(ns_fft); %带噪语音功率
    ns_ps_out(:,n) = ns_ps;
    
%     line(:,1)= ns_ps1(1:160);
%     ns_ps(:,1)=VOICE_FilterLineToBand(line);

%     nn_cut=nn(k:k+Len-1); 
%     nn_win=win.*nn(k:k+Len-1);
%     nn_fft=fft(nn_win,nFFT);
%     nn_out(:,n)=abs(nn_fft);
%     nn_ps(:,n)=nn_fft.*conj(nn_fft); 
    
    if n == 1
      [parametersN,parametersE] = imcra_initialise_parameters(ns_ps); %数据初始化 
    else
      parametersN = estnoise_imcra(ns_ps,parametersN);
      parametersE.noise_ps=parametersN.noise_ps;
      noise_ps_out(:,n) = parametersN.noise_ps;
      Smin_out(:,n) = parametersN.Smin;
      gamma_min(:,n) = ns_ps./(parametersN.Bmin*parametersN.Smin);
%       parametersE.noise_ps= nn_ps;
      parametersE = logmmse_OM(ns_ps,parametersE);
      parametersN.GH1=parametersE.GH1;
    end
    
    noise_ps = parametersN.noise_ps; 
    noise_mu(:,n)=noise_ps;  % magnitude spectrum

    sig =  parametersE.sigout; 
%     line= sig(1:160);
%     mel=VOICE_FilterLineToBand(line);
%     line=VOICE_FilterBandToLin(sig1);
%     sig(1:160,1)=line(1:160);
%     sig(161,1)=0;
%     sig(162:320,1)=flipud(line(2:160,1));
    Gout(:,n)=sig;  % magnitude spectrum
%     gamma_out(:,n)=parametersE.gamma;
%     gammaN_out(:,n)=parametersN.gamma;
%     B=parametersN.Bmin;
%     Bmin_out(:,n)=B;
%     GH1_out0 =  parametersE.GH1; 
%     GH1_out(:,n)=GH1_out0;  % magnitude spectrum
%     smin_sore(:,n) = parametersN.Smin_sw_tild;
    xi_w= ifft(sig .* ns_fft);
    xi_w= real( xi_w);

% --------- Overlap and add ---------------
%
    xfinal(k:k+ninc-1)= x_old+ xi_w(1:Len1);
    x_old= xi_w(Len1+1:Len);  
    k=k+ninc;
end
wavwrite(xfinal,fs,16,'E:\语音增强MATLAB\00_ns_road_ori_matlabsmooth2.wav');
