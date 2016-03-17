echo on

% SPEAKER LOCALIZATION AND TRACKING DEMO
%
% speaker moves from 90?(broadside) to 0?(endfire), back to 90? and
% finally to 180?%
% show usage of program doa_gcc2.m
%
% (c) G. Doblinger, Vienna University of Technology, 2006
% g.doblinger@tuwien.ac.at
% http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%

load mics_harmon.mat            % load microphone array geometry

[X,Fs] = wavread('d:\ÓïÒôÎÄ¼þ\clean\sp01.wav'); % load microphone signals
                                % (8 channel recording, 16 kHz sampling frequency)

d = abs(mics(1,1)-mics(8,1));   % distance of the outmost two microphones

doa_gcc2(1,X(:,1),X(:,8),d);    % SCOT-GCC algorithm, PLEASE WAIT...

pause % PRESS A KEY to continue

% doa_gcc2(2,X(:,1),X(:,8),d);    % PHAT-GCC algorithm, PLEASE WAIT...
