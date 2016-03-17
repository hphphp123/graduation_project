echo on

% ADAPTIVE MICROPHONE ARRAY WITH SPEAKER TRACKING
%
% speaker moves from 90° (broadside) to 0° (endfire), back to 90°, and
% finally to 180°
%
% show usage of program adaptive_beam.m
%
% (c) G. Doblinger, Vienna University of Technology, 2006
% g.doblinger@tuwien.ac.at
% http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%

load mics_harmon.mat            % load microphone array geometry

[X,Fs] = wavread('x_e_16.wav'); % load microphone signals
                                % (8 channel recording, 16 kHz sampling frequency)

y90 = adaptive_beam(X,mics,90); % adaptive array with beamformer look direction = 90°

pause  % PRESS A KEY to continue

y = adaptive_beam(X,mics);      % adaptive array with automatic speaker tracking

% use MATLAB sound playing functions (platform depending)
% to play adaptive array output signals y, and y90
 
