echo on

% SPEAKER LOCALIZATION AND TRACKING DEMO
%
% speaker moves from 90?(broadside) to 0?(endfire), back to 90? and
% finally to 180?%
% show usage of program doa_fastlms.m
%
% (c) G. Doblinger, Vienna University of Technology, 2006
% g.doblinger@tuwien.ac.at
% http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%
% 
% load mics_harmon.mat            % load microphone array geometry
% 
% [X,Fs] = wavread('x_e_16.wav'); % load microphone signals
%                                 % (8 channel recording, 16 kHz sampling frequency)
load data16k.mat
% d = abs(mics(1,1)-mics(8,1));   % distance of the outmost two microphones

doa_fastlms(data(:,1),data(:,2),0.06,512,2048,0.25);    % FASTLMS algorithm, PLEASE WAIT...
