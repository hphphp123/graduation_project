echo on

% SPEAKER LOCALIZATION AND TRACKING DEMO
%
% speaker moves from 90° (broadside) to 0° (endfire), back to 90°, and
% finally to 180°
%
% show usage of program doa_aevd2.m
%
% (c) G. Doblinger, Vienna University of Technology, 2006
% g.doblinger@tuwien.ac.at
% http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%

load mics_harmon.mat            % load microphone array geometry

[X,Fs] = wavread('x_e_16.wav'); % load microphone signals
                                % (8 channel recording, 16 kHz sampling frequency)

d = abs(mics(1,1)-mics(8,1));   % distance of the outmost two microphones

doa_aevd2(X(:,1),X(:,8),d,512,3500,0.3);    % AEVD2 algorithm, PLEASE WAIT...
